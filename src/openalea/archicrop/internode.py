from __future__ import annotations

from math import cos, pi, sin

import numpy as np

import openalea.plantgl.all as pgl


def slim_cylinder(length, radius_base, radius_top):
    """
    Try to construct a cylinder with a low number of triangles.
    """
    rb, rt = radius_base, radius_top
    a1, a2, a3 = 0, 2 * pi / 3.0, 4 * pi / 3.0
    r = rb
    p1 = (r * cos(a1), r * sin(a1), 0)
    p2 = (r * cos(a2), r * sin(a2), 0)
    p3 = (r * cos(a3), r * sin(a3), 0)
    r = rt
    q1 = (r * cos(a1 + pi), r * sin(a1 + pi), length)
    q2 = (r * cos(a2 + pi), r * sin(a2 + pi), length)
    q3 = (r * cos(a3 + pi), r * sin(a3 + pi), length)

    return pgl.TriangleSet(
        [p1, p2, p3, q1, q2, q3],
        [
            (2, 1, 0),
            (3, 4, 5),
            (0, 5, 4),
            (0, 4, 2),
            (2, 4, 3),
            (3, 1, 2),
            (1, 3, 5),
            (5, 0, 1),
        ],
    )


def stem_mesh(length, visible_length, stem_diameter, classic=False, slices=24):  # noqa: ARG001
    """Compute mesh for a stem element
    - classic indicates
    """

    if classic:
        solid = True
        # 6 is the minimal number of slices for a correct computation of star
        #  (percentage error lower than 5)
        slices = 10
        stem = pgl.Tapered(
            stem_diameter / 2.0,
            stem_diameter / 2.0,
            pgl.Cylinder(1.0, visible_length, solid, slices),
        )
        tessel = pgl.Tesselator()
        stem.apply(tessel)
        mesh = tessel.triangulation
    else:
        mesh = slim_cylinder(visible_length, stem_diameter / 2.0, stem_diameter / 2.0)

    return mesh


def stem_dimension(
    h_ins=None,
    d_stem=None,
    internode=None,
    sheath=None,
    d_internode=None,
    d_sheath=None,
    ntop=None,
    plant=1,
):
    """Estimate botanical dimension of stem organs from stem measurements

    Args:
        h_ins: (array) vector of blade insertions height
        d_stem:(float or array) vector of stem diameter
        internode:(array) vector of internode lengths. If None, will be
        estimated using other args
        sheath: (array) vector of sheath lengths. If None, will be estimated
        using other args
        d_internode: (array) vector of intenode diameters. If None, will be
        estimated using other args
        d_sheath: (array) vector of sheath diameters. If None, will be
        estimated using other args
        ntop:(array) vector of leaf position (topmost leaf =1). If None
        (default), stem dimensions are assumed to be from top to base.
        plant: (int or array) vector of plant number

    Returns:
        a dict with estimated sheath and internode dimension
    """

    if h_ins is None and h_ins == internode == sheath:
        h_ins = (60, 50, 40)

    if d_stem is None and d_stem == d_internode == d_sheath:
        d_stem = 0.3

    if h_ins is None:
        sheath = np.array([0] * len(internode)) if sheath is None else np.array(sheath)
        if internode is None:
            internode = np.array([0] * len(sheath))
        else:
            internode = np.array(internode)
        ntop = np.arange(1, len(h_ins) + 1) if ntop is None else np.array(ntop)
        order = np.argsort(-ntop)
        reorder = np.argsort(order)
        h_ins = (internode[order].cumsum() + sheath[order])[reorder]
    else:
        h_ins = np.array(h_ins)
        ntop = np.arange(1, len(h_ins) + 1) if ntop is None else np.array(ntop)
        order = np.argsort(-ntop)
        reorder = np.argsort(order)

        if sheath is None:
            if internode is None:
                sheath = np.array([0] * len(h_ins))
                internode = np.diff([0, *list(h_ins[order])])[reorder]
            else:
                internode = np.array(internode)
                sheath = np.maximum(0, h_ins[order] - internode[order].cumsum())[
                    reorder
                ]

                internode = np.diff([0, *(h_ins[order] - sheath[order]).tolist()])[
                    reorder
                ]
        else:
            sheath = np.array(sheath)
            internode = np.diff([0, *(h_ins[order] - sheath[order]).tolist()])[reorder]

    if d_internode is None:
        d_internode = [d_stem] * len(h_ins) if d_sheath is None else d_sheath
    if d_sheath is None:
        d_sheath = d_internode

    if isinstance(plant, int):
        plant = [plant] * len(ntop)

    return {
            "plant": plant,
            "ntop": ntop,
            "h_ins": h_ins,
            "L_sheath": sheath,
            "W_sheath": d_sheath,
            "L_internode": internode,
            "W_internode": d_internode,
            }


def stem_as_dict(stem_prop, leaf_prop, rank):
    """Return a dictionary with stem properties for a given rank"""
    return {
            "label": "Stem",
            "rank": rank,
            "mature_length": stem_prop["L_internode"][rank-1],
            "length": stem_prop["L_internode"][rank-1],
            "visible_length": stem_prop["L_internode"][rank-1],
            "is_green": True,
            "mature_stem_diameter": stem_prop["W_internode"][rank-1],
            "stem_diameter": stem_prop["W_internode"][rank-1],
            "azimuth": leaf_prop["leaf_azimuth"][rank-1],
            "grow": False,
            "potential_growth_rate": 0.0,
            "age": 0.0,
            "stem_lengths": []
        }