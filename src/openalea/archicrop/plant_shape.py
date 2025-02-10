from __future__ import annotations

import math
import numpy as np
from scipy.integrate import quad
from itertools import product
from scipy.integrate import cumulative_simpson
from scipy.interpolate import splrep

from openalea.mtg.traversal import pre_order2

def geometric_dist(height, nb_phy, q, u0):
    """
    Calculates the heights of leaves'ligules along an axis, according to a geometric series,
    starting from an offset height (u0) = pseudostem height.

    Parameters:
    - height (float): Total height of the plant.
    - nb_phy (int): Number of phytomers (leaves).
    - q (float): Geometric progression factor (controls spacing between leaves).
    - u0 (float): Offset height to start the geometric distribution.

    Returns:
    - List[float]: Normalized distances of leaves along the plant's height.
    """
    if nb_phy <= 0:
        raise ValueError("Number of phytomers (nb_phy) must be greater than zero.")
    if q <= 0:
        raise ValueError("Geometric progression factor (q) must be positive.")
    if u0 >= height:
        raise ValueError("Offset height (u0) must be less than the total height.")

    # Calculate the height available for geometric distribution
    remaining_height = height - u0

    # Generate table of heights for geometric distribution
    table_heights = np.array([i*float(height) / nb_phy if q == 1 else remaining_height * (1 - q**i) / (1 - q**nb_phy) for i in range(1, nb_phy + 1)])
    
    # Add the offset height (u0) to each leaf's position
    normalized_heights = u0 + table_heights
    
    return normalized_heights.tolist()



def collar_heights_kaitaniemi(height, nb_phy):
    """return collar heights"""
    # coefficient k between height of successive collars : 1,17/(1-2,29*exp(-0,79*i))  and k=1,3 for penultimate and k=1,44 for last
    k = [1.17/(1-2.29*math.exp(-0.79*i)) for i in range(1, nb_phy-2)]
    k_pelnultimate = 1.3
    k.append(k_pelnultimate)
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

# def bell_shaped_dist(max_leaf_length, nb_phy, rmax, skew):
#     """returns relative leaf area of individual leaves along bell shaped model"""
#     k = -np.log(skew) * rmax
#     r = np.linspace(1.0 / nb_phy, 1, nb_phy)
#     relative_area = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
#     return relative_area * max_leaf_length

def bell_shaped_dist(Smax, nb_phy, rmax, skew):
    """returns relative leaf area of individual leaves along bell shaped model, so that the sum equals Smax"""
    k = -np.log(skew) * rmax
    r = np.linspace(1.0 / nb_phy, 1, nb_phy)
    relative_area = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    total_area = sum(relative_area)
    normalized_leaf_areas = [area / total_area for area in relative_area]
    return [Smax * la for la in normalized_leaf_areas]




def sigmoid_growth(time, start_tt, end_tt, mature_stem_diameter):
    """
    Sigmoid function for stem diameter growth with inflection at end_tt.

    Parameters:
    - time (float): Current time.
    - start_tt (float): Time when growth starts.
    - end_tt (float): Time when growth reaches the maximum growth rate (inflection point).
    - mature_stem_diameter (float): Maximum diameter.

    Returns:
    - float: Stem diameter at the given time.
    """
    # Adjust growth rate coefficient based on the interval
    k = 4 / (end_tt - start_tt)  # Steeper slope if start and end are close

    # Shift the inflection point to `end_tt`
    t0 = end_tt

    # Compute the sigmoid growth
    diameter = mature_stem_diameter / (1 + math.exp(-k * (time - t0)))
    return diameter



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


def compute_leaf_area_growing_plant(g):
    """returns the leaf area of a growing plant from its MTG(t)"""

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
                    lower_bound = max(L - n.visible_length, 0.0)
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


def compute_height_growing_plant(g):
    """returns the height of a growing plant from its MTG(t)"""

    height = 0

    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))
        for metamer in pre_order2(g, v):
            n = g.node(metamer)
            if n.label is not None:
                if n.label.startswith("Stem"):
                    height += n.visible_length

    return height



def compute_leaf_area_plant_from_params(nb_phy,
                                        Smax,
                                        wl,
                                        rmax,
                                        skew):
    """returns the leaf area of a plant"""

    leaf_areas = []

    leaf_lengths = np.array(bell_shaped_dist(Smax, nb_phy, rmax, skew))

    for L in leaf_lengths:
        blade_area = 2 * wl * (1.51657508881031 - 0.766666666666667) * L**2 # 1.5*wl*L**2 # eg 1.5*0.12*70**2

        leaf_areas.append(blade_area)


    return sum(leaf_areas)






def check_LA_range(params, value_range):
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
        # print(result)
        # Check if the result is within the range
        if min_val <= result <= max_val:
            results.append(values)
    
    return results


def leaf_area(l, L=1, wl=1, alpha=-2.3):
    return 2*l**2*wl*math.sqrt(-alpha) + 2*alpha*l**3*wl/(3*L)
    
def d_leaf_area(l, dl, L=1, wl=1):
    return leaf_area(l+dl,L,wl) - leaf_area(l,L,wl)


def correspondance_dS_dl(wl, L=1):

    dl = 0.0001*L

    l_values = np.arange(0, L, dl)  

    ds_list = []
    for v in l_values:
        ds_list.append(d_leaf_area(v,dl,L,wl))

    corres_dl_dS = {}
    for i,dl in enumerate(l_values):
        corres_dl_dS[dl] = np.cumsum(ds_list)[i]

    return corres_dl_dS


def shape_to_surface(a_leaf, wl):
    _,_,s,r = a_leaf
    r = r[::-1]*wl
    S = cumulative_simpson(x=s, y=r, initial=0)
    S_norm = S/S[-1]
    tck = splrep(x=S_norm, y=s, k=3, task=0)
    return tck


# la_per_plant_stics = max(la_cum)
# print(la_per_plant_stics)
# error = 2000
# value_range = (la_per_plant_stics - error, la_per_plant_stics + error)