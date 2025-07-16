from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import interp1d, splrep
from scipy.optimize import brentq

from .fitting import fit3
from .leaf import leaf_area, leaf_shape_perez


def sr_prevot(nb_segment=100, alpha=-2.3):
    """Leaf shape model from Prevot et al. [year]"""
    beta = -2 * (alpha + np.sqrt(-alpha))
    gamma = 2 * np.sqrt(-alpha) + alpha
    s = np.linspace(0, 1, nb_segment + 1)
    r = alpha * s**2 + beta * s + gamma
    return s, r


def sr_dornbush(nb_segment=100, klig=0.6, swmax=0.55, f1=0.64, f2=0.92):
    """
    Leaf shape model from DornBush et al. 2011
    Parameters
    ----------
    nb_segment: number of points sampled along the leaf to represent the shape
    klig: relative leaf width at leaf base
    swmax: relative distance from leaf tip of the point where leaf width is maximal
    f1 : form factor of the distal leaf segment (near leaf tip)
    f2: form factor if the proximal leaf segment (near leaf base)

    Returns
    -------

    s,r parallel array for curviliear abcissa / relative leaf width along leaf shape
    """

    c1 = 1.0 / f1 - 1
    c2 = brentq(
        lambda x: klig - f2 + (1 - klig) * (1 + 1.0 / x) - 1.0 / np.log(1 + x),
        1,
        1e9,
    )
    ntop = int(nb_segment * swmax)
    offset = 1.0 / nb_segment / 2
    st = np.array([*np.linspace(1 - swmax, 1 - offset, ntop).tolist(), 1])
    rt = np.power((1 - st) / swmax, c1)
    nbase = nb_segment - ntop
    offset = 1.0 / nb_segment / 10
    sb = np.array([0, *np.linspace(offset, 1 - swmax, nbase)[:-1].tolist()])
    rb = klig + (1 - klig) * np.log(1 + c2 * sb / (1 - swmax)) / np.log(1 + c2)
    return np.array(sb.tolist() + st.tolist()), np.array(rb.tolist() + rt.tolist())


def parametric_cereal_leaf(
    nb_segment=10, insertion_angle=50, scurv=0.5, curvature=50, klig=0.6, swmax=0.55, f1=0.64, f2=0.92
):
    """x,y, sr coordinates of points sampling a leaf midrib placed in a vertical plane (origin = leaf base)

    Parameters
    ----------
    insertion_angle: the angle (degree) between stem and leaf at leaf base
    scurv : the relative position on the midrib where 2/3 of total leaf curvature is achieved
    curvature : leaf angular curvature (tip angle - insertion angle, degree)

    Returns
    -------
    x, y coordinates of nb_segments points sampling the leaf midrib
    """
    nseg = min(100, nb_segment)
    x, y = leaf_shape_perez(insertion_angle, scurv, curvature, nseg)
    # s, r = sr_prevot(nseg, alpha)
    s, r = sr_dornbush(nb_segment=nseg, klig=klig, swmax=swmax, f1=f1, f2=f2)
    return fit3(x, y, s, r, nb_points=nb_segment)


def shape_to_surface(leaf, wl):
    _,_,s,r = leaf
    r = r[::-1]*wl
    S = cumulative_simpson(x=s, y=r, initial=0)
    S_norm = S/S[-1]
    return splrep(x=S_norm, y=s, k=3, task=0)


def growing_leaf_area(leaf, length=1, mature_length=1, width_max=1, form_factor=None):
    """
    Leaf area as a function of length
    -------

    Parameters
    ----------
    leaf: x,y,s,r coordinate describing mature leaf shape
    length: current length of the leaf
    mature_length: the length of the leaf once mature
    width_max: maximal width of the mature leaf
    form_factor: (optional) the form_factor of the mature leaf (if known), used to avoid integration

    Returns
    -------
    the area of the leaf corresponding to the distal part up to length (computed with trapeze aera)
    """
    if length >= mature_length and form_factor is not None:
        return length * width_max * form_factor
    _, _, s, r = leaf
    sr_b = 1 - float(length) / mature_length
    return leaf_area(s, r, mature_length, width_max, sr_base=sr_b, sr_top=1)


def truncate_leaf(leaf, fraction=0.1):
    x, y, s, r = leaf
    st = np.linspace(0, fraction, len(s))
    xt = interp1d(s, x)(st)
    yt = interp1d(s, y)(st)
    rt = interp1d(s, r)(1 - fraction + st)
    return xt, yt, st / fraction, rt


def get_base_width(leaf, visible_length=None):
    """Get the width at the base a leaf with a given visibility

    Args:
        leaf:  a x, y, s, r tuple
        length: (float)

    Returns:
        the leaf base width
    """
    _, _, s, r = leaf
    s /= max(s)
    return interp1d(s, r)(1 - visible_length)