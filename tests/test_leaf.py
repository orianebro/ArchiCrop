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

from archicrop.cereals_leaf import (
    leaf_mesh,
    parametric_leaf,
)


def leaf_shape():
    """return x, y, s, r values."""
    return parametric_leaf()


def simple_leaf(leaf, ratio=0.5):
    total_length = 3
    lw_ratio = 8.0

    L_shape = total_length
    Lw_shape = total_length / lw_ratio
    length = total_length * ratio
    s_base = 0
    s_top = 1.0

    return leaf_mesh(leaf, L_shape, Lw_shape, length, s_base, s_top)
