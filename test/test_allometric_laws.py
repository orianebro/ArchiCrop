import matplotlib.pyplot as plt

from openalea.archicrop.plant_shape import bell_shaped_dist

leaf_area = 1
nb_phy = 22
rmax = 0.7
skew = 0.05

leaf_lengths = bell_shaped_dist(leaf_area, nb_phy, rmax, skew)

plt.plot(range(len(leaf_lengths)), leaf_lengths)
plt.show()