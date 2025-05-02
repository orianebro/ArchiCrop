"""Some macro / utilities to run light simpulation on pgl/lpy virtual scene"""

from __future__ import annotations

import numpy as np
import pandas as pd
from alinea.caribu.CaribuScene import CaribuScene

# from alinea.caribu.data_samples import data_path
from alinea.caribu.light import light_sources


def illuminate(scene, light=None, pattern=None, scene_unit="cm", north=0):
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
    if pattern is not None:
        infinite = True
    # if light is not None:
    #     light = light_sources(*light, orientation=north)
    cs = CaribuScene(scene, light=light, scene_unit=scene_unit, pattern=pattern)
    raw, agg = cs.run(direct=True, simplify=True, infinite=infinite)
    return cs, raw["Ei"], raw["Eabs"], pd.DataFrame(agg)


def compute_light_inter(scene, sky, pattern):
    
    # sky = str(data_path('Turtle16soc.light'))
    # zenith = str(data_path('zenith.light'))
    # opts = str(data_path('par_sorghum.opt'))
    # opts = list(map(str, [data_path('par_sorghum.opt'), data_path('nir_sorghum.opt')]))
    materials = {'par': (0.15, 0.06, 0.16, 0.03)}
    # pattern = (-50, -50, 50, 50)
    
    # complete set of files
    cs = CaribuScene(scene=scene, # give mtg rather than scene
                     light=sky, 
                     pattern=pattern, # for toric scene
                     scene_unit='cm',
                     soil_mesh=1, # 1
                     opt=materials) 
    
    # raw,agg=cs.run(simplify=True, infinite=True, sensors={"s1": [(0,0,500),(0,1,500),(1,0,500)]})
    raw,agg=cs.run(simplify=True, infinite=True)

    Qi, Qem, Einc = cs.getIncidentEnergy()
    # print("Qi :", Qi)
    # print("Qem :", Qem)
    # print("Einc :", Einc)
    # Qi is the total horizontal irradiance emitted by sources (m-2)
    # Qem is the sum of the normal-to-the-sources irradiance emitted by sources (m-2)
    # Einc is the total incident energy received on the domain

    # Qi_soil, Einc_soil = cs.getSoilEnergy()
    
    scene_eabs,values_eabs = cs.plot(raw['Eabs'],display=False) 
    scene_ei,values_ei = cs.plot(raw['Ei'],display=False)
    # values = np.array([a / i if i != 0 else 0 for a, i in zip(values_eabs, values_ei)])
    # print("area : ", agg['area'])

    # Compute the 99th percentile 
    v99_eabs = np.percentile(values_eabs, 99)
    v99_ei = np.percentile(values_ei, 99)

    # Convert both lists to numpy arrays
    nvalues_eabs = np.array(values_eabs)
    nvalues_ei = np.array(values_ei)

    # Apply the truncation to both lists
    mask_eabs = nvalues_eabs > v99_eabs  # Identify elements greater than the 99th percentile
    mask_ei = nvalues_ei > v99_ei
    nvalues_eabs[mask_eabs] = v99_eabs
    nvalues_ei[mask_ei] = v99_ei  

    # Calculate the sum of both lists 
    sum_values_eabs = np.sum(nvalues_eabs)
    # sum_values_ei = np.sum(nvalues_ei)  
    # print("Eabs :", sum_values_eabs)
    # print("Ei :", sum_values_ei)


    # Division of sums
    # result = sum_(values_eabs * triangle area) * conv / sum_values_ei if sum_values_ei != 0 else 0
    result = sum_values_eabs / Einc if Einc != 0 else 0 
    # result = (sum_values_ei - Einc_soil) / Einc if Einc != 0 else 0 

    # if result > 1:
    #     PlantGL(scene_eabs, group_by_color=False, property=values_eabs)

    return result, scene_ei, values_ei

# PlantGL(scene, group_by_color=False, property=values)

# Einc is the total incident energy received on the domain
# Eabs (float): the surfacic density of energy absorbed (µmol m⁻² s⁻¹)
# raint (float): from STICS, Photosynthetic Active Radiation intercepted by the canopy 	(MJ m⁻²)
# PAR(MJ m⁻²) = 0.0145 * PAR(µmol m⁻² s⁻¹)

# sum_eabs = 0
# for t in values:
#     sum_eabs += t

# print(sum_eabs*0.0145)