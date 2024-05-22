from __future__ import annotations

import numpy as np
from openalea.mtg.traversal import pre_order2
from scipy.integrate import quad


def geometric_dist(height, nb_phy, q=1):
    """ returns distances between individual leaves along a geometric model """

    if q==1:
        u0=float(height)/nb_phy
    else:
        u0=height*(1.-q)/(1.-q**(nb_phy+1))

    return [u0*q**i for i in range(nb_phy)] # while value < height (but not >> either)


def bell_shaped_dist(max_leaf_length, nb_phy, rmax=.8, skew=0.15):
    """ returns leaf area of individual leaves along bell shaped model """

    k = -np.log(skew) * rmax
    r = np.linspace(1. / nb_phy, 1, nb_phy)
    relative_length = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    # leaf_length = relative_length / relative_length.sum() * max_leaf_length
    leaf_length = relative_length * max_leaf_length
    return leaf_length.tolist()


def compute_leaf_area(g):
    """ returns the leaf area of a plant """

    def scaled_leaf_shape(s, L, alpha=-2.3):
        beta = -2 * (alpha + np.sqrt(-alpha))
        gamma = 2 * np.sqrt(-alpha) + alpha
        r = alpha * (s/L)**2 + beta * (s/L) + gamma
        return r
    
    leaf_areas = []
    # for v in g.vertices():
        # n=g.node(v)
    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))
        for metamer in pre_order2(g, v):
            n = g.node(metamer)
            if n.label is not None:
                if n.label.startswith('Leaf') and n.grow == True:
                    L=n.shape_mature_length
                    alpha=-2.3
                    lower_bound=max(L-n.visible_length, 0.0)
                    upper_bound=L
                    blade_area, error=quad(scaled_leaf_shape, lower_bound, upper_bound, args=(L, alpha))
                    blade_area=2*n.shape_max_width*blade_area
                    leaf_areas.append(blade_area)

                if n.label.startswith('Stem') and n.grow == True:
                    h=n.visible_length
                    radius=n.diameter/2
                    sheath_area=2*np.pi*radius*h
                    leaf_areas.append(sheath_area)

    # filter label
    # g.property('label')

    # PlantGL surface function
                

    return leaf_areas


'''
def sr_prevot(nb_segment=100, alpha=-2.3):
    beta = -2 * (alpha + numpy.sqrt(-alpha))
    gamma = 2 * numpy.sqrt(-alpha) + alpha
    s = numpy.linspace(0, 1, nb_segment + 1)
    r = alpha * s**2 + beta * s + gamma
'''