from __future__ import annotations

from math import pi, sqrt
from operator import itemgetter
from random import random, sample

import numpy as np
from numpy import cos, linspace, sin
from numpy.random import vonmises


def regular(nb_plants, nb_rank, dx, dy, nx=None):
    if nx is None:
        nx = int(nb_plants / nb_rank)
    ny = nb_rank
    domain = ((0, 0), (nx * dx, ny * dy))
    return [
        (i * dx + dx / 2.0, j * dy + dy / 2.0, 0.0)
        for j in range(ny)
        for i in range(nx)
    ], domain


def randomise_position(position, radius):
    az = random() * 2 * np.pi
    r = random() * radius
    dx = r * cos(az)
    dy = r * sin(az)
    x, y, z = position
    return (x + dx, y + dy, z)


def regular_plot(
    inter_plant,
    inter_row,
    nrow,
    plant_per_row,
    noise=0,
    convunit=100,
    center_scene=True,
):
    dx = inter_plant * convunit
    dy = inter_row * convunit
    positions, domain = regular(
        nrow * plant_per_row, int(nrow), dx, dy, int(plant_per_row)
    )
    domain_area = (
        abs(domain[1][0] - domain[0][0])
        / convunit
        * abs(domain[1][1] - domain[0][1])
        / convunit
    )

    # sorting by ranks
    positions = sorted(positions, key=itemgetter(1, 0))
    # add noise
    if noise > 0:
        positions = [randomise_position(x, noise * convunit) for x in positions]
    if center_scene:
        xc = float(domain[1][0] + domain[0][0]) / 2
        yc = float(domain[1][1] + domain[0][1]) / 2
        positions = [(x - xc, y - yc, z) for x, y, z in positions]
        domain = (
            (domain[0][0] - xc, domain[0][1] - yc),
            (domain[1][0] - xc, domain[1][1] - yc),
        )

    return positions, domain, domain_area


def agronomic_plot(
    length=1,
    width=1,
    sowing_density=10,
    inter_row=0.5,
    noise=0,
    convunit=100,
    center_scene=True,
):
    """Returns the number of plants, the positions, the domain (scene units), the domain area (square meter) and the conversion coefficient for meter to scene unit (1/convunit) of a micro-plot specified with agronomical variables
    length (m) is plot dimension along row direction
    width (m) is plot dimension perpendicular to row direction
    sowing density is the density of seeds sawn
    plant_density is the density of plants that are present (after loss due to bad emergence, early death...)
    inter_row (m) is for the  distance between rows
    noise (m) is the radius of the circle where individual plant are randomly positioned
    convunit is the conversion factor from meter to scene unit
    center_scene allows to center the position around origin. If False, the scene is in the x+,y+ sector, the origin being at the lower left corner of the domain

    Rows are parallel to x-axis
    Length and Width are adjusted to produce a canopy centered in its domain and compliant with infinitisation
    """

    # find a (nrow, plant_per_row) sowing design that best match plot dimension
    inter_plant = 1.0 / inter_row / sowing_density
    nrow = max(1, int(round(float(width) / inter_row)))
    plant_per_row = max(1, int(round(float(length) / inter_plant)))
    positions, domain, domain_area = regular_plot(
        inter_plant,
        inter_row,
        nrow,
        plant_per_row,
        noise=noise,
        convunit=convunit,
        center_scene=center_scene,
    )
    n_emerged = int(round(len(positions)))
    positions = sample(positions, n_emerged)

    return n_emerged, positions, domain, domain_area, 1.0 / convunit


def agronomicplotwithdistributions(
    length,
    width,
    sowing_density,
    plant_density,
    inter_row,
    mu=0.0,
    kappa=3.0,
    noise=0,
    convunit=100,
):
    """
    Returns the
        - number of plants
        - the positions
        - azimuthal orientation of each plant
        - the domain
        - the simulated density
        of a micro-plot specified with agronomical variables

    Inputs:
        - length (m) is plot dimension along row direction
        - width (m) is plot dimension perpendicular to row direction
        - sowing density is the density of seeds sawn
        - plant_density is the density of plants that are present (after loss due to bad emergence, early death...)
        - inter_row (m) is for the  distance between rows
        - mu, kappa are the von Mises parameters of the azimuthal distribution
        - noise (%), indicates the precision of the sowing for the inter plant spacing
        - unit (m or cm) is for the unit of the position and domain
    """
    inter_plant = 1.0 / inter_row / sowing_density
    nrow = max(1, int(float(width) / inter_row))
    plant_per_row = max(1, int(float(length) / inter_plant))
    nplants = nrow * plant_per_row
    positions, domain = regular(
        nplants, nrow, inter_plant * convunit, inter_row * convunit
    )
    n_emerged = int(nplants * plant_density / sowing_density)
    positions = sample(positions, n_emerged)
    density = int(
        n_emerged
        / (
            abs(domain[1][0] - domain[0][0])
            / convunit
            * abs(domain[1][1] - domain[0][1])
            / convunit
        )
    )
    azimuths = vonmises(mu, kappa, nplants)
    return n_emerged, positions, azimuths, domain, density


def regularband(nb_plants, nb_rank, dx, dy):
    nx = int(nb_plants / nb_rank)
    ny = nb_rank
    db = sample(range(100), nb_plants)
    domain = ((0, 0), (nx * dx, ny * dy))
    return [
        (
            i * dx + dx / 2.0,
            j * dy + dy / 2.0 + (db[i + j] - 50) / 100.0 * dy / 3.0,
            0.0,
        )
        for j in range(nb_rank)
        for i in range(nx)
    ], domain


def concentric(nb_plants, distance_plant):
    nb_circle = int(sqrt(nb_plants))
    dp = distance_plant
    points = []
    for i in range(nb_circle):
        n = 2 * i + 1
        theta = linspace(0, pi, n)
        radius = i * dp
        x = cos(theta) * radius
        y = sin(theta) * radius
        points.extend([(x[i], y[i], 0.0) for i in range(len(x))])

    return points


def uniform(density, nb_plants):
    return []


def sample_selection(points, gap_fraction):
    """
    Choose a sample from a list of points.
    The gap fraction indicates the proportion of points to remove.
    """
    if gap_fraction == 0:
        return points
    n = len(points)
    k = int(n * (1.0 - gap_fraction / 100.0))
    return sample(points, k)


'''
def clumping_selection(points, nb_clumps, radius_mu, radius_sigma):
    """
    Select a set of points based on andomly chosen groups
    with a radius following a normal distribution.
    """
    n = len(points)
    centers = sample(list(range(n)), nb_clumps)

    for index in centers:
        radius = normalvariate(radius_mu, radius_sigma)
        r2 = radius**2
        # filter to select all the points in a ball of radius r2
        # Add these points to the result
'''


def sample_regular_gaps(points, pattern=[0, 1]):
    """Sample points along pattern.
    Pattern is replicated to become as long as points and then used as a filter (0= gap) on points
    Returns positions of plants and gaps
    """

    if not isinstance(points, list):
        points = [points]

    p = pattern
    length = len(points)
    p = p * (length / len(p)) + [p[i] for i in range(length % len(p))]

    # selection = compress(points,p) only python >= 2.7
    return [point for point, i in zip(points, p) if i], [
        point for point, i in zip(points, p) if not i
    ]


def configuration(Lr, Lx, Ly, scenario_arr, wheat_fraction):
    """
    Return the spatial arrangement of patches.

    Return the spatial arrangement of patches within the field
    and the wheat fraction within each patch.

    Parameters
    ----------
    Lr : int
        number of seasons
    Lx : int
        number of patch along the x-axis
    Ly : int
        number of patch along the y-axis
    scenario_rot : str
        rotation scenario (chose: uniform, random, chessboard, alternate, alternate_rank,
                                alternate_strip2, alternate_strip3, alternate_pairs, etc)
    wheat_fraction : float
        wheat fraction within each patch

    Returns
    -------
    array
    """
    arrangement = np.ones((Lr, Lx, Ly))

    ### MIX OF X% OF WHEAT IN EACH PLOT FOR ALL PLOTS
    if scenario_arr == "uniform":
        arrangement = wheat_fraction * arrangement

    ### RANDOM PATTERN OF X% OF WHEAT
    elif scenario_arr == "random":
        id = np.random.permutation(Lx * Ly)  # noqa: NPY002
        arrangement = arrangement.reshape([Lr * Lx * Ly])
        arrangement[id[0 : int(np.floor(Lx * Ly * (1 - wheat_fraction)))]] = 0
        arrangement = arrangement.reshape([Lr, Lx, Ly])

    ### CHESSBOARD PATTTERN WITH 50% OF WHEAT
    elif scenario_arr == "chessboard":
        for i in range(Lx):
            for j in range(Ly):
                if (i + j) % 2 == 0:
                    arrangement[0, i, j] = 0

    ### ALTERNATE ROWS
    elif scenario_arr == "alternate":
        for i in range(Lx):
            if i % 2 == 1:
                arrangement[:, i, :] = 0

    ### ALTERNATE RANKS (ALONG Ly)
    elif scenario_arr == "alternate_rank":
        for i in range(Ly):
            if i % 2 == 1:
                arrangement[:, :, i] = 0

    ### ALTERNATE STRIPS WITH 2 ROWS PER STRIP
    elif scenario_arr == "alternate_strip2":
        for i in range(Lx):
            if i % 4 == 0:
                arrangement[:, i : i + 2, :] = 0

    ### ALTERNATE STRIPS WITH 3 ROWS PER STRIP
    elif scenario_arr == "alternate_strip3":
        for i in range(Lx):
            if i % 6 == 0:
                arrangement[:, i : i + 3, :] = 0

    ### ALTERNATE PAIRS OF PATCH ALONG Lx
    elif scenario_arr == "alternate_pairs":
        for i in range(Lx):
            for j in range(Ly):
                if i % 2 == 0 and j % 4 == 0:
                    arrangement[:, i, j : j + 2] = 0
                    arrangement[:, i + 1, j + 2 : j + 3] = 0

    ### ALTERNATE DOUBLE-PAIRS ALONG Lx
    elif scenario_arr == "alternate_doublepairs":
        for i in range(Lx):
            if i % 2 == 1:
                arrangement[:, i, :] = 0
        for i in range(Lx):
            for j in range(Ly):
                if i % 2 == 0 and j % 8 == 0:
                    arrangement[:, i, j : j + 4] = 0
                    arrangement[:, i + 1, j : j + 4] = 1

    ### TWO SUB-FIELDS
    elif scenario_arr == "two_subs":
        for i in range(Lx):
            if i < (Lx / 2):
                arrangement[:, i, :] = 0

    ### FOUR SUB-FIELDS
    elif scenario_arr == "four_subs":
        for i in range(Lx):
            for j in range(Ly):
                if i < (Lx / 2) and j < (Ly / 2):
                    arrangement[:, i, j] = 0
                if i >= (Lx / 2) and j >= (Ly / 2):
                    arrangement[:, i, j] = 0

    ### 4-PATCH SQUARE
    elif scenario_arr == "4_square":
        for i in range(Lx):
            for j in range(Ly):
                if i % 4 == 0 and j % 4 == 0:
                    arrangement[:, i : i + 2, j : j + 2] = 0
                    arrangement[:, i + 2 : i + 4, j + 2 : j + 4] = 0

    ### 9-PATCH SQUARES
    elif scenario_arr == "9_square":
        for i in range(Lx):
            for j in range(Ly):
                if i % 6 == 0 and j % 6 == 0:
                    arrangement[:, i : i + 3, j : j + 3] = 0
                    arrangement[:, i + 3 : i + 6, j + 3 : j + 6] = 0

    ############ TEST DIFFERENT PROPORTIONS ACCORDING TO ARRANGEMENT ############
    ### 1/5 IN ALTERNATE ROWS
    elif scenario_arr == "1_5_alternate":
        for i in range(Lx):
            if i % 5 == 0:
                arrangement[:, i, :] = 0

    ### 2/5 IN ALTERNATE ROWS
    elif scenario_arr == "2_5_alternate":
        for i in range(Lx):
            for j in range(Lx):
                if i % 5 == 0:
                    arrangement[:, i, :] = 0
                if j % 5 == 2:
                    arrangement[:, j, :] = 0

    ### 3/5 IN ALTERNATE ROWS
    elif scenario_arr == "3_5_alternate":
        for i in range(Lx):
            for j in range(Lx):
                if i % 5 == 0:
                    arrangement[:, i, :] = 0
                if j % 5 == 2:
                    arrangement[:, j : j + 2, :] = 0

    ### 4/5 IN ALTERNATE ROWS
    elif scenario_arr == "4_5_alternate":
        for i in range(Lx):
            if i % 5 == 1:
                arrangement[:, i : i + 4, :] = 0

    ### 1/5 IN 4-PATCH SQUARE
    elif scenario_arr == "1_5_square":
        for i in range(Lx):
            for j in range(Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i : i + 2, j : j + 2] = 0
                    arrangement[:, i + 2 : i + 4, j + 4 : j + 6] = 0
                    arrangement[:, i + 4 : i + 6, j + 8 : j + 10] = 0
                    arrangement[:, i + 6 : i + 8, j + 2 : j + 4] = 0
                    arrangement[:, i + 8 : i + 10, j + 6 : j + 8] = 0

    ### 2/5 IN 4-PATCH SQUARE
    elif scenario_arr == "2_5_square":
        for i in range(Lx):
            for j in range(Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i : i + 2, j : j + 2] = 0
                    arrangement[:, i : i + 2, j + 4 : j + 6] = 0
                    arrangement[:, i + 2 : i + 4, j + 2 : j + 4] = 0
                    arrangement[:, i + 2 : i + 4, j + 6 : j + 8] = 0
                    arrangement[:, i + 4 : i + 6, j + 4 : j + 6] = 0
                    arrangement[:, i + 4 : i + 6, j + 8 : j + 10] = 0
                    arrangement[:, i + 6 : i + 8, j : j + 2] = 0
                    arrangement[:, i + 6 : i + 8, j + 6 : j + 8] = 0
                    arrangement[:, i + 8 : i + 10, j + 2 : j + 4] = 0
                    arrangement[:, i + 8 : i + 10, j + 8 : j + 10] = 0

    ### 3/5 IN 4-PATCH SQUARE
    elif scenario_arr == "3_5_square":
        for i in range(Lx):
            for j in range(Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i : i + 2, j : j + 2] = 0
                    arrangement[:, i : i + 2, j + 4 : j + 8] = 0
                    arrangement[:, i + 2 : i + 4, j + 2 : j + 4] = 0
                    arrangement[:, i + 2 : i + 4, j + 6 : j + 10] = 0
                    arrangement[:, i + 4 : i + 6, j : j + 2] = 0
                    arrangement[:, i + 4 : i + 6, j + 4 : j + 6] = 0
                    arrangement[:, i + 4 : i + 6, j + 8 : j + 10] = 0
                    arrangement[:, i + 6 : i + 8, j : j + 4] = 0
                    arrangement[:, i + 6 : i + 8, j + 6 : j + 8] = 0
                    arrangement[:, i + 8 : i + 10, j + 2 : j + 6] = 0
                    arrangement[:, i + 8 : i + 10, j + 8 : j + 10] = 0

    ### 4/5 IN 4-PATCH SQUARE
    elif scenario_arr == "4_5_square":
        for i in range(Lx):
            for j in range(Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i : i + 2, j : j + 8] = 0
                    arrangement[:, i + 2 : i + 4, j + 2 : j + 10] = 0
                    arrangement[:, i + 4 : i + 6, j + 4 : j + 10] = 0
                    arrangement[:, i + 4 : i + 6, j : j + 2] = 0
                    arrangement[:, i + 6 : i + 8, j : j + 4] = 0
                    arrangement[:, i + 6 : i + 8, j + 6 : j + 10] = 0
                    arrangement[:, i + 8 : i + 10, j : j + 6] = 0
                    arrangement[:, i + 8 : i + 10, j + 8 : j + 10] = 0

    return arrangement
