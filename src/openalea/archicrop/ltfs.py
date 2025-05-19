"""Some macro / utilities to run light simpulation on pgl/lpy virtual scene """
import pandas
import numpy
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.light import light_sources  # here to avoid import line in notebook
'''
from openalea.lpy import Lsystem


def run_lsystem(lsystem='leafy.lpy', parameters=None):
    if parameters:
        lsys = Lsystem(lsystem, {'parameters': parameters})
    else:
        lsys = Lsystem(lsystem)
    lstring = lsys.iterate()
    lscene = lsys.sceneInterpretation(lstring)
    return lsys, lstring, lscene
'''

def illuminate(scene, light=None, pattern=None, scene_unit='cm', north=0, labels=None):
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
    df = pandas.DataFrame(agg)
    if labels is not None:
        labs = pandas.DataFrame(labels)
        df = labs.merge(df.reset_index(), left_index=True, right_index=True)
    return cs, raw['Ei'], df


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
            raise ValueError('At least one grid dimension (dx, dy) should be specified')
    return (-0.5 * dx, -0.5 * dy,
             0.5 * dx, 0.5 * dy)
