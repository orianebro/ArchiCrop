"""Some macro / utilities to run light simpulation on pgl/lpy virtual scene """
import pandas
import numpy
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.light import light_sources  # here to avoid import line in notebook


def illuminate(scene, light=None, pattern=None, scene_unit='cm', north=0):
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
    cs = CaribuScene(scene, light=light,scene_unit=scene_unit, pattern=pattern)
    raw, agg = cs.run(direct=True, simplify=True, infinite=infinite)
    return cs, raw['Ei'], pandas.DataFrame(agg)


