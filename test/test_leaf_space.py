from __future__ import annotations

from openalea.archicrop.cereals_leaf import parametric_leaf
from openalea.archicrop.fitting import mesh4, plantgl_shape
from openalea.archicrop.geometry import arrange_leaf
from openalea.plantgl.all import *


def make_leaf(inclination):
    leaf = arrange_leaf(leaf=parametric_leaf(), inclination=inclination)
    return plantgl_shape(*mesh4(
        leaf, length_max=10, length=10, s_base=0, s_top=1, radius_max=1, volume=0
    )) 


