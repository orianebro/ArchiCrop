from __future__ import annotations

import math
import numpy as np
from scipy.integrate import quad
from itertools import product

from openalea.mtg.traversal import pre_order2


def geometric_dist(height, nb_phy, q=1, u0=None):
    """returns distances between individual leaves along a geometric model"""

    if u0 is None:
        if q == 1:
            u0 = float(height) / nb_phy
        else:
            u0 = height * (1.0 - q) / (1.0 - q ** (nb_phy + 1))

    # print(u0*(1-q**(nb_phy))/(1-q))

    table_heights = np.array([i * u0 * q**i for i in range(1,nb_phy+1)])


    return [i * height / table_heights[-1] for i in table_heights]

def collar_heights_kaitaniemi(height, nb_phy):
    """return collar heights"""
    # coefficient k between height of successive collars : 1,17/(1-2,29*exp(-0,79*i))  and k=1,3 for penultimate and k=1,44 for last
    k = [1.17/(1-2.29*math.exp(-0.79*i)) for i in range(1, nb_phy-2)]
    k_pelnultinate = 1.3
    k.append(k_pelnultinate)
    k_flag = 1.44
    k.append(k_flag)

    # Compute first collar height
    k_product = 1
    k_sum = 1
    for j in k:
        k_product *= j
        k_sum += k_product

    u1 = height * 1/k_sum

    # Compute other collars height
    internode_lengths = [u1]
    k_product = 1  
    for f in k:
        k_product *= f
        internode_lengths.append(k_product * u1)

    collar_heights = np.cumsum(np.array(internode_lengths))

    return collar_heights


def bell_shaped_dist(max_leaf_length, nb_phy, rmax=0.8, skew=0.15):
    """returns leaf area of individual leaves along bell shaped model"""

    k = -np.log(skew) * rmax
    r = np.linspace(1.0 / nb_phy, 1, nb_phy)
    relative_length = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    # leaf_length = relative_length / relative_length.sum() * max_leaf_length
    leaf_length = relative_length * max_leaf_length
    return leaf_length.tolist()


def compute_leaf_area(g):
    """returns the leaf area of a plant"""

    def scaled_leaf_shape(s, L, alpha=-2.3):
        beta = -2 * (alpha + np.sqrt(-alpha))
        gamma = 2 * np.sqrt(-alpha) + alpha
        r = alpha * (s / L) ** 2 + beta * (s / L) + gamma
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
                if n.label.startswith("Leaf") and n.grow == True:
                    L = n.mature_length
                    wl = 0.12
                    blade_area = 2 * (1.51657508881031*L**2*wl - 0.766666666666667*L**2*wl)
                    # alpha = -2.3
                    # lower_bound = max(L - n.visible_length, 0.0)
                    # upper_bound = L
                    # blade_area, error = quad(
                    #     scaled_leaf_shape, lower_bound, upper_bound, args=(L, alpha)
                    # )
                    # blade_area = 2 * n.shape_max_width * blade_area
                    leaf_areas.append(blade_area)

                # if n.label.startswith("Stem") and n.grow == True:
                #     h = n.visible_length
                #     radius = n.diameter / 2
                #     sheath_area = 2 * np.pi * radius * h
                #     leaf_areas.append(sheath_area)

    # filter label
    # g.property('label')

    # PlantGL surface function

    return leaf_areas


def compute_leaf_area_pot_plant(g):
    """returns the leaf area of a potential plant"""

    def scaled_leaf_shape(s, L, alpha=-2.3):
        beta = -2 * (alpha + np.sqrt(-alpha))
        gamma = 2 * np.sqrt(-alpha) + alpha
        r = alpha * (s / L) ** 2 + beta * (s / L) + gamma
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
                if n.label.startswith("Leaf"):
                    L = n.mature_length
                    alpha = -2.3
                    lower_bound = max(L - n.mature_length, 0.0)
                    upper_bound = L
                    blade_area, error = quad(
                        scaled_leaf_shape, lower_bound, upper_bound, args=(L, alpha)
                    )
                    blade_area = 2 * n.shape_max_width * blade_area
                    leaf_areas.append(blade_area)

                # if n.label.startswith("Stem") and n.grow == True:
                #     h = n.visible_length
                #     radius = n.diameter / 2
                #     sheath_area = 2 * np.pi * radius * h
                #     leaf_areas.append(sheath_area)

    # filter label
    # g.property('label')

    # PlantGL surface function

    return leaf_areas


"""
def sr_prevot(nb_segment=100, alpha=-2.3):
    beta = -2 * (alpha + numpy.sqrt(-alpha))
    gamma = 2 * numpy.sqrt(-alpha) + alpha
    s = numpy.linspace(0, 1, nb_segment + 1)
    r = alpha * s**2 + beta * s + gamma
"""

def compute_leaf_area_plant_from_params(nb_phy,
                                        max_leaf_length,
                                        wl,
                                        rmax,
                                        skew):
    """returns the leaf area of a plant"""

    leaf_areas = []

    leaf_lengths = np.array(bell_shaped_dist(max_leaf_length=max_leaf_length, nb_phy=nb_phy, rmax=rmax, skew=skew))

    for L in leaf_lengths:
        blade_area = 2 * wl * (1.51657508881031 - 0.766666666666667) * L**2 # 1.5*wl*L**2 # eg 1.5*0.12*70**2

        leaf_areas.append(blade_area)


    return sum(leaf_areas)




def check_la_range(params, value_range):
    """
    Computes function values for all parameter combinations
    and keeps only those whose results fall within a given range.

    :param params: A dictionary where keys are parameter names
                   and values are lists of possible values.
                   Example: {'x': [1, 2], 'y': [3, 4], 'z': [5, 6]}
    :param value_range: A tuple (min_val, max_val) defining the range.
    :return: A list of tuples (combination, function_value).
    """
    min_val, max_val = value_range

    # Generate all possible combinations of parameters
    param_names = list(params.keys())
    param_names_for_la = ["nb_phy", "max_leaf_length", "wl", "rmax", "skew"]
    # param_values = list(params.values())
    param_values = [
        v if isinstance(v, list) else [v]
        for v in params.values()
    ]
    param_values_for_la = [param for i,param in enumerate(param_values) if param_names[i] in param_names_for_la]

    
    combinations = list(product(*param_values_for_la))
    
    # Filter combinations based on the range
    results = []
    for combination in combinations:
        # Convert the combination into a dictionary
        values = dict(zip(param_names_for_la, combination))
        # Calculate the function value
        result = compute_leaf_area_plant_from_params(**values)
        print(result)
        # Check if the result is within the range
        if min_val <= result <= max_val:
            results.append((values, result))
    
    return results

# la_per_plant_stics = max(la_cum)
# print(la_per_plant_stics)
# error = 2000
# value_range = (la_per_plant_stics - error, la_per_plant_stics + error)