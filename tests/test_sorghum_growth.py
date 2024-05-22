from __future__ import annotations

import numpy as np

tmax = 100


def plant_height(t):
    return 300 / (1 + np.exp(-0.1 * (t - tmax / 2)))


internodes = []
internodes.append(0)
phyllochron = 5
i = 0

for t in range(tmax):
    internodes[i] += plant_height(t) - plant_height(
        t - 1
    )  # or internodes[i] = plant_height(t)
    if t % phyllochron == 0:
        internodes.append(internodes[i])
        i += 1


print(internodes)
