import numpy as np
import matplotlib.pyplot as plt

# from openalea.archicrop.plant_shape import geometric_dist

def geometric_dist(height, nb_phy, q, u0):
    """
    Calculates the distances between individual leaves along a geometric progression model,
    starting from an offset height (u0).

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



def build_stem(nb_phy = 20, height = 222, stem_q = 1.1, nb_young_phy = 5, young_phy_height = 2):

    pseudostem_height = nb_young_phy * young_phy_height

    pseudostem = np.array([young_phy_height*i for i in range(1, nb_young_phy+1)])
    stem = np.array(geometric_dist(height, nb_phy-nb_young_phy, q=stem_q, u0=pseudostem_height))
    insertion_heights = np.concatenate((pseudostem, stem), axis=0)

    for i in range(1,len(insertion_heights)):
        assert(insertion_heights[i-1] < insertion_heights[i])

    return insertion_heights


nb_phy = 20

for stem_q in np.arange(1.0,1.4,0.1):
    print(stem_q)
    insertion_heights = build_stem(nb_phy = 20, height = 222, stem_q = stem_q, nb_young_phy = 5, young_phy_height = 2)
    plt.plot(range(1,nb_phy+1), insertion_heights)

plt.show()