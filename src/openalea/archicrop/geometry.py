from __future__ import annotations

import openalea.plantgl.all as pgl

from .internode import stem_mesh
from .leaf import leaf_mesh


def mesh_area(mesh):
    if mesh:
        sc = pgl.SurfComputer(pgl.Discretizer())
        mesh.apply(sc)
        return sc.surface
    return None


def _is_iterable(x):
    try:
        iter(x)
    except TypeError:
        return False
    return True


def as_tuples(pgl_3List, offset=0):
    """return pgl list of 3 numbers kind (indes3, vector3) as a list of python
    tuples
    """
    if not _is_iterable(offset):
        offset = [offset] * 3
    return [(i[0] + offset[0], i[1] + offset[1], i[2] + offset[2]) for i in pgl_3List]


def addSets(pglset1, pglset2, translate=(0, 0, 0)):
    """
    Create a new TriangleSet by addition of two existing ones
    if translate is not None, pglset2 is translated with vector translate
    """
    points = as_tuples(pglset1.pointList) + as_tuples(
        pglset2.pointList, offset=translate
    )
    index = as_tuples(pglset1.indexList) + as_tuples(
        pglset2.indexList, offset=len(pglset1.pointList)
    )
    return pgl.TriangleSet(points, index)


def compute_organ_geometry(element_node, classic=False):
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""
    n = element_node
    geom = None

    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
            if n.shape is not None and n.srb is not None:
                geom = leaf_mesh(
                    n.shape,
                    n.mature_length,
                    n.shape_max_width,
                    n.visible_length,
                    n.srb,
                    n.srt,
                    # flipx allows x-> -x to place the shape along
                    #  with the tiller positioned with
                    # turtle.down()
                    flipx=True,
                    stem_diameter=n.stem_diameter,
                )
            if n.lrolled > 0:
                rolled = stem_mesh(
                    n.lrolled, n.lrolled, n.d_rolled, classic
                )
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        geom = stem_mesh(n.length, n.visible_length, n.stem_diameter, n.stem_diameter, classic)

    return geom


### for growing organs

def compute_growing_organ_geometry(element_node, classic=False):
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""

    n = element_node
    geom = None

    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
            if n.shape is not None and n.srb is not None:
                # stem_diameter=0
                geom = leaf_mesh(
                    leaf=n.shape,
                    L_shape=n.mature_length,
                    Lw_shape=n.shape_max_width,
                    length=n.visible_length,
                    s_base=n.srb,
                    s_top=n.srt,
                    # flipx allows x-> -x to place the shape along
                    # with the tiller positioned with
                    # turtle.down()
                    flipx=True,
                    inclination=n.inclination, 
                    stem_diameter=n.stem_diameter,
                )

            ### Add a mesh that is senescent
            if n.senescent_length > 0.0001: # filter less than 0.001 mm leaves  # noqa: SIM102
                if n.srt is not None: # and n.shape is not None
                    geom_senescent = leaf_mesh(
                        leaf=n.shape,
                        L_shape=n.mature_length,
                        Lw_shape=n.shape_max_width,
                        length=n.visible_length,
                        s_base=n.srt,
                        s_top=1.,
                        # flipx allows x-> -x to place the shape along
                        # with the tiller positioned with
                        # turtle.down()
                        flipx=True,
                        inclination=n.inclination, 
                        stem_diameter=n.stem_diameter,
                    )
                    n.geometry_senescent = geom_senescent


            if n.lrolled > 0:
                rolled = stem_mesh(
                    n.lrolled, n.lrolled, n.d_rolled, classic
                )
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        geom = stem_mesh(length=n.length, visible_length=n.visible_length, stem_diameter=n.stem_diameter, classic=classic)

    return geom



