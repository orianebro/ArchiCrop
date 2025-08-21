from __future__ import annotations

from math import cos, pi, radians, sin

import numpy as np
from scipy.integrate import simpson, trapezoid

import openalea.plantgl.all as pgl

from .fitting import mesh4, plantgl_shape

# def leaf_shape_perez(nb_segment=100, insertion_angle=50, delta_angle=180, coef_curv=-0.2):
#     def _curvature(s, coef_curv):
#         return ((1 + coef_curv) * (s**2)) / (1 + coef_curv * (s**2))
#         # positionRelativeLeaf=vector of relative position on the leaf [0;1]
#         # decli_ini = declination (=inclination from the vertical) at leaf insertion (degree) [0;180]
#         # decli_final = declination (=inclination from the vertical) at leaf tip (degree) [decli_ini; 180]
#         # coefCurv= coefficient of leaf curvature [-1,inf] -1=no curvature inf=curvature at insertion
#         # leafLength = length of the leaf
#
#     s = numpy.linspace(0,1,nb_segment+1)
#     ds = 1. / (nb_segment)
#     angle_simu = _curvature(s, coef_curv=coef_curv) * numpy.radians(
#         delta_angle) + numpy.radians(insertion_angle)
#     dx = numpy.array([0] + (ds * numpy.sin(angle_simu)).tolist())[:-1]
#     dy = numpy.array([0] + (ds * numpy.cos(angle_simu)).tolist())[:-1]
#     x, y = numpy.cumsum(dx), numpy.cumsum(dy)
#     length = numpy.sum(numpy.sqrt(dx**2 + dy**2))
#     return x / length, y / length


def leaf_shape_perez(insertion_angle=50, scurv=0.5, curvature=50, nb_segment=100):
    """x,y coordinates of points sampling a leaf midrib placed in a vertical plane (origin = leaf base)

    Parameters
    ----------
    insertion_angle: the angle (degree) between stem and leaf at leaf base
    scurv : the relative position on the midrib where 2/3 of total leaf curvature is achieved
    curvature : leaf angular curvature (tip angle - insertion angle, degree)

    Returns
    -------
    x, y coordinates of nb_segments points sampling the leaf midrib
    """

    def _curvature(s, coef_curv):
        return ((1 + coef_curv) * (s**2)) / (1 + coef_curv * (s**2))

    # inputs

    s = np.linspace(0, 1, nb_segment + 1)
    # ds = 1.0 / (nb_segment)
    # fraction of delta angle reach at l
    frac_l = 2.0 / 3
    # curvature coefficient of the first section (before l)
    coefCurv_1 = -0.2
    # curvature coefficient of the second section (after l)
    coefCurv_2 = 5

    # leaf tip angle
    tip_angle = insertion_angle + curvature

    # leaf angle at l
    l_angle = insertion_angle + frac_l * (tip_angle - insertion_angle)

    # angles in the first section of the curve
    angle_simu_1 = _curvature(s, coef_curv=coefCurv_1) * np.radians(
        l_angle - insertion_angle
    ) + np.radians(insertion_angle)

    # angles in the second section of the curve
    angle_simu_2 = _curvature(s[1:], coef_curv=coefCurv_2) * np.radians(
        tip_angle - l_angle
    ) + np.radians(l_angle)

    # all angles
    angle_simu = np.array(angle_simu_1.tolist() + angle_simu_2.tolist())
    coef_l = [scurv] * len(s) + [1 - scurv] * len(s[1:])

    dx = np.array([0, *(coef_l * np.sin(angle_simu)).tolist()])[:-1]
    dy = np.array([0, *(coef_l * np.cos(angle_simu)).tolist()])[:-1]
    x, y = np.cumsum(dx), np.cumsum(dy)
    length = np.sum(np.sqrt(dx**2 + dy**2))
    return x / length, y / length


def leaf_morpho_rel(nb_segment=10, w0=0.5, lm=0.5):
    a0 = w0
    c0 = (w0 - 1) / (lm**2)
    b0 = -2 * c0 * lm

    c1 = -1 / (1 - lm) ** 2
    b1 = -2 * c1 * lm
    a1 = -b1 - c1

    s = np.linspace(0, 1, nb_segment + 1)

    r1 = np.array(a0 + b0 * s[s <= lm] + c0 * s[s <= lm] ** 2)
    r2 = np.array(a1 + b1 * s[s > lm] + c1 * s[s > lm] ** 2)
    r = np.concatenate([r1, r2])
    return s, r


def leaf_area(s, r, Lshape=1, Lwshape=1, sr_base=0, sr_top=1):
    """surface of a blade element, positioned with two relative curvilinear absisca"""

    sr_base = min([1, max([0, sr_base])])
    sr_top = min([1, max([sr_base, sr_top])])
    sre = [sr for sr in zip(s, r) if (sr_base < sr[0] < sr_top)]
    if len(sre) > 0:
        se, re = zip(*sre)
        snew = [sr_base, *list(se), sr_top]
        rnew = [np.interp(sr_base, s, r), *list(re), np.interp(sr_top, s, r)]
    else:
        snew = [sr_base, sr_top]
        rnew = [np.interp(sr_base, s, r), np.interp(sr_top, s, r)]

    S = trapezoid(rnew, snew) * Lshape * Lwshape

    return S  # noqa: RET504



def form_factor(leaf, simps=False):
    """Compute the form factor of a leaf"""
    _, _, s, r = leaf
    if simps:
        return simpson(r, s)    
    return leaf_area(s, r, 1, 1, 0, 1)



def arrange_leaf(leaf, stem_diameter=0, inclination=1, relative=True):
    """Arrange a leaf to be placed along a stem with a given inclination.

    Args:
        leaf: a x, y, s, r tuple describing leaf shape
        stem_diameter: the diameter of the sem at the leaf insertion point
        inclination: if relative=False, the leaf basal inclination (deg). A
        multiplier to leaf basal inclination angle otherwise
        relative: (bool) controls the meaning of inclination parameter

    Returns:
        a modified x, y, s, r tuple

    """

    x, y, s, r = list(map(np.array, leaf))
    if relative and inclination == 1:
        x1, y1 = x, y
    else:
        basal_inclination = pgl.angle((x[1] - x[0], y[1] - y[0]), (0, 1))

        if relative:
            angle = inclination * basal_inclination
            angle = min(pi, angle)
        else:
            angle = radians(inclination)

        rotation_angle = basal_inclination - angle

        # rotation of the midrib
        cos_a = cos(rotation_angle)
        sin_a = sin(rotation_angle)

        x1 = x[0] + cos_a * x - sin_a * y
        y1 = y[0] + sin_a * x + cos_a * y

    leaf = x1 + stem_diameter / 2.0, y1, s, r

    return leaf  # noqa: RET504


def leaf_mesh(
    leaf,
    L_shape=1,
    Lw_shape=1,
    length=1,
    s_base=0,
    s_top=1,
    flipx=False,
    twist=0,
    volume=0,
    stem_diameter=0,
    inclination=1,
    relative=True,
):
    """Compute mesh for a leaf element along a scaled leaf shape

    Args:
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

    Returns:
        a PlanGl mesh representing the element
    """
    shape = arrange_leaf(
        leaf,
        stem_diameter=stem_diameter / L_shape,
        inclination=inclination,
        relative=relative,
    )

    # flip to position leaves along tiller emitted
    if flipx:
        # to position leaves along tiller emitted
        shape = (-shape[0],) + shape[1:]  # noqa: RUF005

    mesh = mesh4(
        shape, L_shape, length, s_base, s_top, Lw_shape, twist=twist, volume=volume
    )

    if mesh:
        pts, ind = mesh
        mesh = None if len(ind) < 1 else plantgl_shape(pts, ind)
    else:
        if length > 0:
            msg = f"Error : no mesh. s_base = {s_base}, s_top = {s_top}, length = {length}"
            raise ValueError(msg)
        mesh = None

    return mesh


def leaf_dimension(
    area=None, length=None, width=None, ntop=None, form_factor=None, plant=1, wl=0.2
):
    """Estimate blade dimension and/or compatibility with leaf shapes form
    factors

    Args:
        area: (array) vector of blade area. If None, will be estimated using
         other args
        length: (array) vector of blade lengths. If None, will be estimated
        using other args
        width: (array) vector of blade widths. If None, will be estimated using
         other args
        ntop: (array) vector of leaf position (topmost leaf =1). If None
        (default), leaf dimensions are assumed to be from top to base.
        form_factor: (object) The (width * length) / Surface ratio If None
         (default) the adel default shape will be used
        plant: (int or array) vector of plant number
        wl: (float) the width / length ratio used to estimates dimensions in
        case of incomplete data

    Returns:
        a dict with estimated blade dimensions

    """

    if area is None and length is None and width is None:
        area = (15, 20, 30)

    form_factor = np.array(0.75) if form_factor is None else np.array(form_factor)

    wl = np.array(wl)

    if area is None:
        if length is None:
            width = np.array(width)
            length = width / np.array(wl)
        elif width is None:
            length = np.array(length)
            width = length * np.array(wl)
        else:
            length = np.array(length)
            width = np.array(width)
        ntop = np.arange(1, len(length) + 1) if ntop is None else np.array(ntop)
        area = form_factor * length * width
    else:
        area = np.array(area)
        ntop = np.arange(1, len(area) + 1) if ntop is None else np.array(ntop)
        # adjust length/width if one is  None or overwrite width if all are set
        if length is None:
            if width is None:
                length = np.sqrt(area / form_factor / wl) # do with scaled max curvilinear abscissa ?
                width = length * wl
            else:
                width = np.array(width)
                length = area / form_factor / width
        else:
            length = np.array(length)
            width = area / form_factor / length

    if isinstance(plant, int):
        plant = [plant] * len(ntop)

    return {
            "plant": plant,
            "ntop": ntop,
            "wl": wl,
            "L_blade": length,
            "W_blade": width,
            "S_blade": area,
            "form_factor": form_factor,
            }



def leaf_as_dict(stem_prop, leaf_prop, rank, wl):
    """Return a dictionary with leaf properties for a given rank"""
    return {
            "label": "Leaf",
            "rank": rank,
            "shape": leaf_prop["leaf_shape"][rank-1],
            "mature_length": leaf_prop["L_blade"][rank-1],
            "length": leaf_prop["L_blade"][rank-1],
            "visible_length": leaf_prop["L_blade"][rank-1],
            "leaf_area": leaf_prop["S_blade"][rank-1],
            "visible_leaf_area": leaf_prop["S_blade"][rank-1],
            "senescent_area": 0.0,
            "senescent_length": 0.0,
            # "form_factor": leaf_prop["form_factor"][rank-1],
            "wl": wl,
            "tck": leaf_prop["tck"][rank-1], 
            "is_green": True,
            "srb": 0,
            "srt": 1,
            "lrolled": 0,
            "d_rolled": 0,
            "shape_max_width": leaf_prop["W_blade"][rank-1],
            "mature_stem_diameter": stem_prop["W_internode"][rank-1],
            "stem_diameter": stem_prop["W_internode"][rank-1],
            "grow": False,
            "potential_growth_rate": 0.0,
            "dead": False,
            "age": 0.0,
            "inclination": 0.0,
            "leaf_lengths": [0.0],
            "senescent_lengths": [0.0]
        }

