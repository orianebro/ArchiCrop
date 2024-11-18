"""Some macro / utilities to run light simpulation on pgl/lpy virtual scene"""

from __future__ import annotations

import numpy as np
import pandas as pd
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.data_samples import data_path
from alinea.caribu.light import light_sources  # here to avoid import line in notebook


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
    if light is not None:
        light = light_sources(*light, orientation=north)
    cs = CaribuScene(scene, light=light, scene_unit=scene_unit, pattern=pattern)
    raw, agg = cs.run(direct=True, simplify=True, infinite=infinite)
    return cs, raw["Ei"], pd.DataFrame(agg)


def compute_light_inter(scene):
    
    sky = str(data_path('Turtle16soc.light'))
    # zenith = str(data_path('zenith.light'))
    # opts = map(str, [data_path('par.opt'), data_path('nir.opt')])
    pattern = (-50, -50, 50, 50)
    
    # complete set of files
    cs = CaribuScene(scene=scene, # give mtg rather than scene
                     light=sky, 
                     pattern=pattern, # for toric scene
                     scene_unit='cm',
                     soil_mesh=1) 
    
    raw,agg=cs.run(simplify=True)
    
    scene,values = cs.plot(raw['Eabs'],display=False)
    
    v99 = np.percentile(values, 99)
    nvalues=np.array(values)
    nvalues[nvalues>v99]=v99
    values = nvalues.tolist()
    return values

# PlantGL(scene, group_by_color=False, property=values)

# Eabs (float): the surfacic density of energy absorbed (µmol m⁻² s⁻¹)
# raint (float): from STICS, Photosynthetic Active Radiation intercepted by the canopy 	(MJ m−2)
# PAR(MJ m-2) = 0.0145 * PAR(µmol m⁻² s⁻¹)

# sum_eabs = 0
# for t in values:
#     sum_eabs += t

# print(sum_eabs*0.0145)