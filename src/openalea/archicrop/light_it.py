"""Some macro / utilities to run light simpulation on pgl/lpy virtual scene"""

from __future__ import annotations

import numpy as np
import pandas as pd
from alinea.caribu.CaribuScene import CaribuScene

from openalea.archicrop.display import build_scene
from openalea.archicrop.stand import compute_domain
from openalea.archicrop.stics_io import stics_weather_3d
from openalea.astk.sky_irradiance import sky_irradiance
from openalea.astk.sky_sources import caribu_light_sources, sky_sources

# Einc is the total incident energy received on the domain
# Eabs (float): the surfacic density of energy absorbed (MJ m⁻² s⁻¹)
# raint (float): from STICS, Photosynthetic Active Radiation intercepted by the canopy 	(MJ m⁻²) --> faPAR !!!

def illuminate(scene, light=None, domain=None, scene_unit='cm', labels=None, direct=True): # north=0):
    """Illuminate scene

    Args:
        scene: the scene (plantgl)
        light: lights. If None a vertical light is used
        pattern: the toric canopy pattern. If None, no pattern is used
        scene_unit: string indicating length unit in the scene (eg 'cm')
        north: the angle (deg, positive clockwise) from X+ to
         North (default: 0)

    Returns:

    """
    infinite = False
    if domain is not None:
        infinite = True
    # if light is not None:
    #     light = light_sources(*light, orientation=north)
    cs = CaribuScene(scene, light=light,scene_unit=scene_unit, pattern=domain, soil_mesh=1)
    raw, agg = cs.run(simplify=True, infinite=infinite, direct=direct) # direct=True
    df = pd.DataFrame(agg)  # noqa: PD901
    if labels is not None:
        labs = pd.DataFrame(labels)
        df = labs.merge(df.reset_index(), left_index=True, right_index=True)  # noqa: PD901
    return cs, raw['Eabs'], df


def mean_leaf_irradiance(df):
    df['Energy'] = df['Eabs'] * df['area']
    agg = df.loc[(df.is_green) & (df.label=='Leaf'),('plant','Energy','area')].groupby('plant').agg('sum')
    agg['Irradiance'] = agg['Energy'] / agg['area']
    return agg


def toric_canopy_pattern(dx=80, dy=5, density=None):
    if density is not None:
        if dx is not None:
            dy = 1. / density / (dx / 100.) * 100
        elif dy is not None:
            dx = 1. / density / (dy / 100.) * 100
        else:
            raise ValueError('At least one grid dimension (dx, dy) should be specified')  # noqa: EM101
    return (-0.5 * dx, -0.5 * dy,
             0.5 * dx, 0.5 * dy)


def compute_light_inter(scene, sky, pattern):
    
    # sky = str(data_path('Turtle16soc.light'))
    # zenith = str(data_path('zenith.light'))
    # opts = str(data_path('par_sorghum.opt'))
    # opts = list(map(str, [data_path('par_sorghum.opt'), data_path('nir_sorghum.opt')]))
    materials = {'par': (0.15, 0.06, 0.16, 0.03)}
    
    # complete set of files
    cs = CaribuScene(scene=scene, # give mtg rather than scene
                     light=sky, 
                     pattern=pattern, # for toric scene
                     scene_unit='cm',
                     soil_mesh=1, # 1
                     opt=materials) 
    
    raw,agg=cs.run(simplify=True, infinite=True)

    Qi, Qem, Einc = cs.getIncidentEnergy()
    # Qi is the total horizontal irradiance emitted by sources (m-2)
    # Qem is the sum of the normal-to-the-sources irradiance emitted by sources (m-2)
    # Einc is the total incident energy received on the domain

    # Qi_soil, Einc_soil = cs.getSoilEnergy()
    
    scene_eabs,values_eabs = cs.plot(raw['Eabs'],display=False) 

    v99_eabs = np.percentile(values_eabs, 99)
    nvalues_eabs = np.array(values_eabs)
    mask_eabs = nvalues_eabs > v99_eabs  
    nvalues_eabs[mask_eabs] = v99_eabs 
    sum_values_eabs = np.sum(nvalues_eabs)

    # Division of sums
    result = sum_values_eabs / Einc if Einc != 0 else 0 
    # result = (sum_values_ei - Einc_soil) / Einc if Einc != 0 else 0 

    return result  # noqa: RET504


# Viewer.display(scene_eabs, group_by_color=False, property=values_eabs)
# PlantGL(scene, group_by_color=False, property=values)


location = {  
    'longitude': 3.87,
    'latitude': 45,
    'altitude': 56,
    'timezone': 'Europe/Paris'}

def light_interception(weather_file, daily_dynamics, sowing_density, location, mtgs, zenith=False, save_scenes=False, inter_row=70):
    '''Compute light interception on plants with fitting parameters
    Args:
        weather_file: path to the weather file
        daily_dynamics: daily dynamics from STICS
        sowing_density: sowing density
        location: dictionary with location parameters (longitude, latitude, altitude, timezone)
        mtgs: list of MTGs for each plant
        zenith: if True, use zenith light sources
        save_scenes: if True, save the scenes as images
    Returns:
        nrj_per_plant: list of energy per plant
    '''
    # Read weather data
    df_weather = stics_weather_3d(filename=weather_file, daily_dynamics=daily_dynamics)

    # Sowing pattern
    domain = compute_domain(sowing_density, inter_row = inter_row) # cm

    # Define incident PAR
    if zenith:
        par_incident = [value["Incident PAR"] for value in daily_dynamics.values()]
    else:
        par_incident = list(df_weather.itertuples())

    # Compute light interception for each plant at each time step
    # par_caribu = []
    # nrj_per_leaf = []
    nrj_per_plant = {}
    # For each plant
    for k, mtgs_plant in mtgs.items():
        if mtgs_plant[0] is None:
            nrj_per_plant[k] = [None] * len(par_incident)
        else:
            # aggs_tmp = []
            nrj_tmp = []
            # For each time step
            for i,(mtg, par) in enumerate(zip(mtgs_plant, par_incident)):
                # if i%1==0:
                # Compute light sources
                if zenith:
                    lights = [(par,(0,0,-1))]
                else:
                    irr = sky_irradiance(daydate=par.daydate, day_ghi=par.rad, **location)
                    sun, sky = sky_sources(sky_type='clear_sky', sky_irradiance=irr, scale='global')
                    lights = caribu_light_sources(sun, sky)
                # Build and illuminate scene
                scene, labels = build_scene(mtg, (0,0,0), senescence=False)
                cs, raw, agg = illuminate(scene, light=lights, labels=labels, domain=domain, direct=False) # --> cf PARaggregators in caribu scene node
                # Compute energy per leaf
                df_mod = mean_leaf_irradiance(agg)  # noqa: F841
                nrj_tmp.append(agg.loc[agg['label'] == 'Leaf']['Energy'].values)
                # aggs_tmp.append(agg)
                # Save scene if required
                if save_scenes:
                    scene_tmp = cs.plot(raw, display=False)[0]
                    scene_tmp.save(f'scene_{i}.png') # not as images !!!

            nrj_per_plant[k] = [sum(growing_plant) for growing_plant in nrj_tmp]
            # par_caribu.append(aggs_tmp)

    '''
    # Calculate energy per leaf and irradiance per plant
    nrj_per_leaf = []
    # irr_per_plant = []

    for illuminated_plant in par_caribu:
        nrj_tmp = []
        # irr_tmp = []
        for df_scene in illuminated_plant:
            df_mod = mean_leaf_irradiance(df_scene)  # noqa: F841
            nrj_tmp.append(df_scene.loc[df_scene['label'] == 'Leaf']['Energy'].values)
            # irr_tmp.append(df_mod['Irradiance'].values[0])  
        nrj_per_leaf.append(nrj_tmp)
        # irr_per_plant.append(irr_tmp)

    nrj_per_plant = [[sum(growing_plant) for growing_plant in plant] for plant in nrj_per_leaf] # to dict !!!!
    '''

    return nrj_per_plant 
