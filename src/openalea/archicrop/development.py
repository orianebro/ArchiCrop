from __future__ import annotations


def add_development(g, vid, tt, dtt, rate):
    """
    Add dynamic properties on the mtg to simulate development

    :param g: MTG, MTG of a plant
    :param vid: int, vertex id of the plant component
    :param tt: float, current thermal time (in °C.day)
    :param dtt: float, thermal time increment (in °C.day)
    :param rate: float, rate of development (in °C.day)

    :return: tt
    """
    g.node(vid).start_tt = tt
    g.node(vid).end_tt = tt + dtt
    tt += rate
    return tt


def compute_potential_growth_rate(g, vid):
    """
    Compute the potential growth rate in length of an organ based on its elongation time.
    :param g: MTG, MTG of a plant
    :param vid: int, vertex id of the organ
    """
    if g.node(vid).label.startswith("Leaf"):
        g.node(vid).potential_growth_rate = g.node(vid).leaf_area / (g.node(vid).end_tt - g.node(vid).start_tt)
    elif g.node(vid).label.startswith("Stem"):
        g.node(vid).potential_growth_rate = g.node(vid).mature_length / (g.node(vid).end_tt - g.node(vid).start_tt)
