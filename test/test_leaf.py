"""Test leaf shape.

Build parameters for a leaf
1. Build a leaf shape (x, y, s, r)
2. Set the params
        leaf: a x,y,s,r tuple describing leaf shape
        L_shape: length of the shape
        Lw_shape: width of the shape
        length: the total visible length to be meshed
        s_base: normalised position (on the shape) of the start of the element
        s_top: normalised position (on the shape) of the end of the element
        flipx: should leaf shape be flipped ?
        twist:
        volume: float value of the thickness of the leaf.
              Default is 0. Else it indicates the relative depth of the leaf
              along the midrib.
        stem_diameter: the diameter of the sem at the leaf insertion point
        inclination: if relative=False, the leaf basal inclination (deg). A
        multiplier to leaf basal inclination angle otherwise
        relative: (bool) controls the meaning of inclination parameter

3. Understand how to build a growing leaf for 0 to 1;
"""

from __future__ import annotations

from pytest import approx  # noqa: PT013

from openalea.archicrop.cereal_leaf import parametric_cereal_leaf
from openalea.archicrop.leaf import leaf_mesh
from openalea.plantgl.all import surface


def leaf_shape():
    """return x, y, s, r values."""
    return parametric_cereal_leaf()


def simple_leaf(leaf, ratio=1):
    total_length = 3
    lw_ratio = 8.0

    L_shape = total_length
    Lw_shape = total_length / lw_ratio
    length = total_length * ratio
    s_base = 0
    s_top = 1.0

    return leaf_mesh(leaf, L_shape, Lw_shape, length, s_base, s_top)

def test_leaf_surface():
    sh = leaf_shape()
    mesh1 = simple_leaf(sh, ratio=0.001)
    assert(surface(mesh1) == approx(0., abs=1e-4))

    total_surface = surface(simple_leaf(sh,ratio=1))

    mesh2 = simple_leaf(sh, ratio=1/3)
    surf2 = surface(mesh2)

    assert(surf2 == approx(0.1566, abs=1e-4))
    assert(0. <= surf2 <= total_surface)

    mesh3 = simple_leaf(sh, ratio=2/3)
    surf3 = surface(mesh3)

    assert(surf3 == approx(0.5024, abs=1e-4))
    assert(0. <= surf3 <= total_surface)

test_leaf_surface()