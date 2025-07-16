from __future__ import annotations

from itertools import cycle

import numpy as np


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
        + (np.random.random(size) - 0.5) * 2 * phyllotactic_deviation  # noqa: NPY002
    )
    azim = azim % 360
    return np.where(azim <= 180, azim, azim - 360)


