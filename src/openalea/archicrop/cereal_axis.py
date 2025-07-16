from __future__ import annotations

import math

import numpy as np


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
        msg = "Number of phytomers (nb_phy) must be greater than zero."
        raise ValueError(msg)
    if q <= 0:
        msg = "Geometric progression factor (q) must be positive."
        raise ValueError(msg)
    if u0 >= height:
        msg = "Offset height (u0) must be less than the total height."
        raise ValueError(msg)

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

    return np.cumsum(np.array(internode_lengths))


# def bell_shaped_dist(max_leaf_length, nb_phy, rmax, skew):
#     """returns relative leaf area of individual leaves along bell shaped model"""
#     k = -np.log(skew) * rmax
#     r = np.linspace(1.0 / nb_phy, 1, nb_phy)
#     relative_area = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
#     return relative_area * max_leaf_length

def bell_shaped_dist(leaf_area, nb_phy, rmax, skew):
    """returns relative leaf area of individual leaves along bell shaped model, so that the sum equals leaf_area"""
    k = -np.log(skew) * rmax
    r = np.linspace(1.0 / nb_phy, 1, nb_phy)
    relative_area = np.exp(-k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    total_area = sum(relative_area)
    normalized_leaf_areas = [area / total_area for area in relative_area]
    return [leaf_area * la for la in normalized_leaf_areas]



