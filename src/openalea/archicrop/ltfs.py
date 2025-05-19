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


def leaf_irradiance(df, lstring, leaf_name='Leaf', aggregate=False):
    leaves = [pid for pid,m in enumerate(lstring) if m.name == leaf_name] # adapt this to MTG
    if len(leaves) > 0:
        df = df.loc[df.index.isin(leaves), :]
        df.loc[:, 'leaf_rank'] = range(1, len(leaves) + 1)
    if aggregate:
        df = df.apply(numpy.mean)
    return df


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
