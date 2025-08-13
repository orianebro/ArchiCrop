from __future__ import annotations

# import matplotlib.pyplot as plt
import numpy as np

from openalea.archicrop.cereal_axis import geometric_dist


def build_stem(nb_phy = 20, height = 222, stem_q = 1.1, nb_young_phy = 5, young_phy_height = 2):

    pseudostem_height = nb_young_phy * young_phy_height

    pseudostem = np.array([young_phy_height*i for i in range(1, nb_young_phy+1)])
    stem = np.array(geometric_dist(height, nb_phy-nb_young_phy, q=stem_q, u0=pseudostem_height))
    insertion_heights = np.concatenate((pseudostem, stem), axis=0)

    for i in range(1,len(insertion_heights)):
        assert(insertion_heights[i-1] < insertion_heights[i])

    # return insertion_heights

'''
nb_phy = 20

for stem_q in np.arange(1.0,1.4,0.1):
    print(stem_q)
    insertion_heights = build_stem(nb_phy = 20, height = 222, stem_q = stem_q, nb_young_phy = 5, young_phy_height = 2)
    plt.plot(range(1,nb_phy+1), insertion_heights)

plt.show()
'''