"""Helpers for designing cereals plants"""

from __future__ import annotations

import base64
from itertools import cycle

import numpy as np
import pandas as pd
from scipy.integrate import simpson
from scipy.interpolate import interp1d


def json_numpy_obj_hook(dct):
    """
    Decodes a previously encoded numpy ndarray
    with proper shape and dtype
    :param dct: (dict) json encoded ndarray
    :return: (ndarray) if input was an encoded ndarray
    """
    if isinstance(dct, dict) and "__ndarray__" in dct:
        data = base64.b64decode(dct["__ndarray__"])
        return np.frombuffer(data, dct["dtype"]).reshape(dct["shape"])
    return dct


def get_form_factor(leaf):
    """
    return form factor for a x,y,s,r tuple referencing a leaf
    """
    _, _, s, r = leaf
    return simpson(r, s)


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


def blade_dimension(
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
        a pandas dataframe with estimated blade dimensions

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

    return pd.DataFrame(
        {
            "plant": plant,
            "ntop": ntop,
            "wl": wl,
            "L_blade": length,
            "W_blade": width,
            "S_blade": area,
            "form_factor": form_factor,
        }
    )


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
        a pandas dataframe with estimated sheath and internode dimension
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

    return pd.DataFrame(
        {
            "plant": plant,
            "ntop": ntop,
            "h_ins": h_ins,
            "L_sheath": sheath,
            "W_sheath": d_sheath,
            "L_internode": internode,
            "W_internode": d_internode,
        }
    )


def leaf_azimuth(
    size=1,
    phyllotactic_angle=180,
    phyllotactic_deviation=15,
    plant_orientation=0,
    spiral=False,
):
    """Generate leaf azimuth series

    Args:
        size: the size of the sample
        phyllotactic_angle: if spiral=False (default) the phyllotactic angle (deg) bet
        ween 'left and right' leaves. If spiral is True, the angle between
        successive leaves (deg)
        phyllotactic_deviation: half-amplitude of deviation around phyllotactic
        angle (deg)
        plant_orientation : first azimuth of the series (deg, from X+ positive
        counter-clockwise)

    Returns:
        an array of azimuths (deg, from X+, positive counter-clockwise)
    """
    if size == 1:
        return plant_orientation
    if spiral:
        main = np.arange(0, size) * phyllotactic_angle
    else:
        it = cycle((0, phyllotactic_angle))
        main = np.array([next(it) for i in range(size)])
    azim = (
        plant_orientation
        + main
        + (np.random.random(size) - 0.5) * 2 * phyllotactic_deviation
    )
    azim = azim % 360
    return np.where(azim <= 180, azim, azim - 360)
