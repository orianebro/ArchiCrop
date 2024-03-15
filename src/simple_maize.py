""" A collection of simple allometric functions for parameterising the
architecture of maize"""
import numpy
from itertools import cycle
from scipy.interpolate import interp1d

from plant_design import get_form_factor, blade_dimension, \
    stem_dimension, load_leaf_db
from fitting import fit3

def maize_leaves_path():
    import os
    import openalea.cereals
    data_dir = os.path.dirname(openalea.cereals.__file__)
    return os.path.join(data_dir, 'leaves_simple.db')


def bell_shaped_dist(plant_area=1, nb_phy=15, rmax=.7, skew=0.15):
    """ returns leaf area of individual leaves along bell shaped model """

    k = -numpy.log(skew) * rmax
    r = numpy.linspace(1. / nb_phy, 1, nb_phy)
    relative_surface = numpy.exp(
        -k / rmax * (2 * (r - rmax) ** 2 + (r - rmax) ** 3))
    leaf_area = relative_surface / relative_surface.sum() * plant_area
    return leaf_area.tolist()


def geometric_dist(height=15, nb_phy=15, q=1):
    """ returns distances between individual leaves along a geometric model """

    if q == 1:
        u0 = float(height) / nb_phy
    else:
        u0 = height * (1. - q) / (1. - q ** (nb_phy + 1))

    return [u0 * q ** i for i in range(nb_phy)]


def leaf_azimuth(size=1, phyllotactic_angle=180, phyllotactic_deviation=15, plant_orientation=0, spiral=False):
    """ Generate leaf azimuth series

    Args:
        size: the size of the sample
        phyllotactic_angle: if spiral=False (default) the phyllotactic angle (deg) bet
        ween 'left and right' leaves. If spiral is True, the angle between
        successive leaves (deg)
        phyllotactic_deviation: half-amplitude of deviation around phyllotactic
        angle (deg)
        plant_orientation : first azimuth of the serie (deg, from X+ positive
        counter-clockwise)

    Returns:
        an array of azimuths (deg, from X+, positive counter-clockwise)
    """
    if size == 1:
        return plant_orientation
    if spiral:
        main = numpy.arange(0, size) * phyllotactic_angle
    else:
        it = cycle((0, phyllotactic_angle))
        main = numpy.array([next(it) for i in range(size)])
    azim = plant_orientation + main + (numpy.random.random(size) - 0.5) * 2 * phyllotactic_deviation
    azim = azim % 360
    return numpy.where(azim <= 180, azim, azim - 360)


def leaf_shape_perez(nb_segment=100, insertion_angle=50, delta_angle=180, coef_curv=-0.2):
    def _curvature(s, coef_curv):
        return ((1 + coef_curv) * (s**2)) / (1 + coef_curv * (s**2))
        # positionRelativeLeaf=vector of relative position on the leaf [0;1]
        # decli_ini = declination (=inclination from the vertical) at leaf insertion (degree) [0;180]
        # decli_final = declination (=inclination from the vertical) at leaf tip (degree) [decli_ini; 180]
        # coefCurv= coefficient of leaf curvature [-1,inf] -1=no curvature inf=curvature at insertion
        # leafLength = length of the leaf

    s = numpy.linspace(0,1,nb_segment+1)
    ds = 1. / (nb_segment)
    angle_simu = _curvature(s, coef_curv=coef_curv) * numpy.radians(
        delta_angle) + numpy.radians(insertion_angle)
    dx = numpy.array([0] + (ds * numpy.sin(angle_simu)).tolist())[:-1]
    dy = numpy.array([0] + (ds * numpy.cos(angle_simu)).tolist())[:-1]
    x, y = numpy.cumsum(dx), numpy.cumsum(dy)
    length = numpy.sum(numpy.sqrt(dx**2 + dy**2))
    return x / length, y / length


def sr_prevot(nb_segment=10, alpha=-2.5):
    beta = -2 * (alpha + numpy.sqrt(-alpha))
    gamma = 2 * numpy.sqrt(-alpha) + alpha
    s = numpy.linspace(0, 1, nb_segment + 1)
    r = alpha * s**2 + beta * s + gamma
    return s, r


def parametric_leaf(nb_segment=10, insertion_angle=50, delta_angle=180,
                    coef_curv=-0.2, alpha=-2.5):
    """

    Args:
        nb_segment:
        insertion_angle:
        delta_angle:
        coef_curv:
        alpha:

    Returns:

    """
    nseg = min(100, nb_segment)
    x, y = leaf_shape_perez(nseg, insertion_angle, delta_angle, coef_curv)
    s, r = sr_prevot(nseg, alpha)
    return fit3(x, y, s, r, nb_points=nb_segment)


def simple_maize(plant_area=10000, plant_height=200, pseudostem_height=20,
                 phytomer=16, rmax=0.67, pseudostem_dist=1.4, stem_dist=1.,
                 diam_base=2.5, diam_top=1, leaves=None, phyllotactic_angle=180,
                 phyllotactic_deviation=15, plant_orientation=0, wl=0.1, skew=0.15, seed=None):
    """ Generate a detailed parameter set for maize simulation from global
    parameters

    Args:
        plant_area: the total leafy area of the plant
        plant_height:  height of the stem
        pseudostem_height:  height of the highestpseudo-stem collar
        phytomer: number of phytomers
        leaf rank function
        rmax: control the magitude of the leaf area as a function of
        leaf rank function
        pseudostem_dist: parameter contrling spacing between pseudostem leaves
        stem_dist: parameter controling spacing between stem leaves
        diam_base: diameter at the base of the stem
        diam_top: dimaeter at the top of the stem
        leaves: a {rank: (x, y, s, r)} dict of tuples defining the (x,y)
        coordinates of leaf midribs and (s,r) coordinate defining the radius as
        a function of distance to leaf base
        phyllotactic_angle: phyllotactic angle between successive leaves (deg)
        phyllotactic_deviation: absolute deviation around phyllotactic angle (deg)
        plant_orientation: azimuth orientation (deg, from X+, positive
        counter-clockwise) of the first leaf of the plant

    Returns: A pandas Dataframe with individual organ dimensions and geometric
    parameters needed by cereals constructors

    """

    numpy.random.seed(seed)
    ranks = range(1, phytomer + 1)
    ntop = max(ranks) - numpy.array(ranks) + 1
    if leaves is None:
        path = maize_leaves_path()
        db = load_leaf_db(path)
        leaves = {rank: db.get(str(rank), db['10'])[0] for rank in ranks}

    phytomer = int(phytomer)

    # Lejeune an Bernier formula + col =
    nb_young_phy = int(round((phytomer - 1.95) / 1.84 / 1.3))

    # compute the leaf surface
    leaf_area = numpy.array(
        bell_shaped_dist(plant_area=plant_area, nb_phy=phytomer, rmax=rmax,
                         skew=skew))

    # distances between leaves
    pseudostem = geometric_dist(pseudostem_height, nb_young_phy,
                                pseudostem_dist)
    stem = geometric_dist(plant_height - pseudostem_height,
                          phytomer - nb_young_phy, stem_dist)
    internode = pseudostem + stem
    # stem diameters
    diameter = ([diam_base] * nb_young_phy +
                numpy.linspace(diam_base, diam_top,
                               phytomer - nb_young_phy).tolist())

    ff = [get_form_factor(leaves[rank]) for rank in ranks]
    blades = blade_dimension(area=leaf_area, form_factor=ff, ntop=ntop, wl=wl)
    stem = stem_dimension(internode=internode, d_internode=diameter, ntop=ntop)
    df = blades.merge(stem)

    df['leaf_azimuth'] = leaf_azimuth(size=len(ranks), phyllotactic_angle=phyllotactic_angle, phyllotactic_deviation=phyllotactic_deviation,
                                      plant_orientation=plant_orientation)
    df['leaf_rank'] = ranks
    df['leaf_shape'] = [leaves[n] for n in df.leaf_rank]

    return df


