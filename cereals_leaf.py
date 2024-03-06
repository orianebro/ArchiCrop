""" Parametric leaf used for simple maize"""
import numpy
from math import pi, cos, sin, radians
import fitting as fitting
import openalea.plantgl.all as pgl
from scipy.integrate import trapz
from scipy.optimize import brentq



# def leaf_shape_perez(nb_segment=100, insertion_angle=50, delta_angle=180, coef_curv=-0.2):
#     def _curvature(s, coef_curv):
#         return ((1 + coef_curv) * (s**2)) / (1 + coef_curv * (s**2))
#         # positionRelativeLeaf=vector of relative position on the leaf [0;1]
#         # decli_ini = declination (=inclination from the vertical) at leaf insertion (degree) [0;180]
#         # decli_final = declination (=inclination from the vertical) at leaf tip (degree) [decli_ini; 180]
#         # coefCurv= coefficient of leaf curvature [-1,inf] -1=no curvature inf=curvature at insertion
#         # leafLength = length of the leaf
#
#     s = numpy.linspace(0,1,nb_segment+1)
#     ds = 1. / (nb_segment)
#     angle_simu = _curvature(s, coef_curv=coef_curv) * numpy.radians(
#         delta_angle) + numpy.radians(insertion_angle)
#     dx = numpy.array([0] + (ds * numpy.sin(angle_simu)).tolist())[:-1]
#     dy = numpy.array([0] + (ds * numpy.cos(angle_simu)).tolist())[:-1]
#     x, y = numpy.cumsum(dx), numpy.cumsum(dy)
#     length = numpy.sum(numpy.sqrt(dx**2 + dy**2))
#     return x / length, y / length



def leaf_shape_perez(insertion_angle=50, scurv=0.5, curvature=50, nb_segment = 100):
    """ x,y coordinates of points sampling a leaf midrib placed in a vertical plane (origin = leaf base)

    Parameters
    ----------
    insertion_angle: the angle (degree) between stem and leaf at leaf base
    scurv : the relative position on the midrib where 2/3 of total leaf curvature is achieved
    curvature : leaf angular curvature (tip angle - insertion angle, degree)

    Returns
    -------
    x, y coordinates of nb_segments points sampling the leaf midrib
    """

    def _curvature(s, coef_curv):
        return ((1 + coef_curv) * (s**2)) / (1 + coef_curv * (s**2))

    #inputs

    s = numpy.linspace(0,1,nb_segment+1)
    ds = 1. / (nb_segment)
    # fraction of delta angle reach at l
    frac_l = 2. / 3
    # curvature coefficient of the first section (before l)
    coefCurv_1 = -0.2
    # curvature coefficient of the second section (after l)
    coefCurv_2 = 5

    #leaf tip angle
    tip_angle = insertion_angle + curvature

    # leaf angle at l
    l_angle = insertion_angle + frac_l*(tip_angle - insertion_angle)

    #angles in the first section of the curve
    angle_simu_1 = _curvature(s, coef_curv=coefCurv_1) * numpy.radians(l_angle-insertion_angle) + numpy.radians(insertion_angle)

    # angles in the second section of the curve
    angle_simu_2 = _curvature(s[1:], coef_curv=coefCurv_2) * numpy.radians(tip_angle - l_angle) + numpy.radians(l_angle)

    #all angles
    angle_simu=numpy.array(angle_simu_1.tolist()+angle_simu_2.tolist())
    coef_l=[scurv]*len(s)+[1-scurv]*len(s[1:])

    dx = numpy.array([0] + (coef_l * numpy.sin(angle_simu)).tolist())[:-1]
    dy = numpy.array([0] + (coef_l * numpy.cos(angle_simu)).tolist())[:-1]
    x, y = numpy.cumsum(dx), numpy.cumsum(dy)
    length = numpy.sum(numpy.sqrt(dx**2 + dy**2))
    return x / length, y / length


def sr_prevot(nb_segment=100, alpha=-2.3):
    beta = -2 * (alpha + numpy.sqrt(-alpha))
    gamma = 2 * numpy.sqrt(-alpha) + alpha
    s = numpy.linspace(0, 1, nb_segment + 1)
    r = alpha * s**2 + beta * s + gamma
    return s, r


def sr_dornbush(nb_segment=100, klig=0.6, swmax=0.55, f1=0.64, f2=0.92):
    """
    Leaf shape model from DornBush et al. 2011
    Parameters
    ----------
    nb_segment: number of points sampled along the leaf to represent the shape
    klig: relative leaf width at leaf base
    swmax: relative distance from leaf tip of the point where leaf width is maximal
    f1 : form factor of the distal leaf segment (near leaf tip)
    f2: form factor if the proximal leaf segment (near leaf base)

    Returns
    -------

    s,r parralel array for curviliear abcissa / relative leaf width along leaf shape
    """

    c1 = 1. / f1 - 1
    c2 = brentq(lambda x: klig - f2 + (1 -klig) * (1 + 1. / x) - 1. / numpy.log(1 + x), 1, 1e9)
    ntop = int(nb_segment * swmax)
    offset = 1. / nb_segment / 2
    st = numpy.array(numpy.linspace(1 - swmax, 1 - offset, ntop).tolist() + [1])
    rt = numpy.power((1 - st) / swmax, c1)
    nbase = nb_segment - ntop
    offset = 1. / nb_segment / 10
    sb = numpy.array([0] + numpy.linspace(offset, 1 - swmax, nbase)[:-1].tolist())
    rb = klig + (1 - klig) * numpy.log(1 + c2 * sb / (1 - swmax)) / numpy.log(1 + c2)
    return numpy.array(sb.tolist() + st.tolist()), numpy.array(rb.tolist() + rt.tolist())


def leaf_morpho_rel(nb_segment=10, w0=0.5, lm=0.5):
    a0 = w0
    c0 = (w0 - 1) / (lm ** 2)
    b0 = -2 * c0 * lm

    c1 = -1 / (1 - lm) ** 2
    b1 = -2 * c1 * lm
    a1 = -b1 - c1

    s = numpy.linspace(0, 1, nb_segment + 1)

    r1 = numpy.array(a0+b0*s[s <=lm]+c0*s[s <=lm]**2)
    r2 = numpy.array(a1+b1*s[s >lm]+c1*s[s >lm]**2)
    r = numpy.concatenate([r1,r2])
    return s, r

# def leaf_area_rel(r,w0,lm):
#     a0 = w0
#     c0 = (w0 - 1) / (lm ** 2)
#     b0 = -2 * c0 * lm
#
#     c1 = -1 / (1 - lm) ** 2
#     b1 = -2 * c1 * lm
#     a1 = -b1 - c1
#
#     if (r<=lm): s = (a0 * r + b0 / 2 * r ** 2 + c0 / 3 * r ** 3)
#     if (r > lm):s= a0 * lm + b0 / 2 * lm ** 2 + c0 / 3 * lm ** 3 + (a1 * r + b1 / 2 * r ** 2 + c1 / 3 * r ** 3) - (a1 * lm + b1 / 2 * lm ** 2 + c1 / 3 * lm ** 3)
#
#     return s, r


def parametric_leaf(nb_segment=10, insertion_angle=50, scurv=0.5,
                     curvature=50, alpha=-2.3):
    """ x,y, sr coordinates of points sampling a leaf midrib placed in a vertical plane (origin = leaf base)

    Parameters
    ----------
    insertion_angle: the angle (degree) between stem and leaf at leaf base
    scurv : the relative position on the midrib where 2/3 of total leaf curvature is achieved
    curvature : leaf angular curvature (tip angle - insertion angle, degree)

    Returns
    -------
    x, y coordinates of nb_segments points sampling the leaf midrib
    """
    nseg = min(100, nb_segment)
    x, y = leaf_shape_perez(insertion_angle, scurv, curvature, nseg)
    s, r = sr_prevot(nseg, alpha)
    return fitting.fit3(x, y, s, r, nb_points=nb_segment)

def blade_elt_area(s, r, Lshape=1, Lwshape=1, sr_base=0, sr_top=1):
    """ surface of a blade element, positioned with two relative curvilinear absisca"""

    S = 0
    sr_base = min([1, max([0, sr_base])])
    sr_top = min([1, max([sr_base, sr_top])])
    sre = [sr for sr in zip(s, r) if (sr_base < sr[0] < sr_top)]
    if len(sre) > 0:
        se, re = zip(*sre)
        snew = [sr_base] + list(se) + [sr_top]
        rnew = [numpy.interp(sr_base, s, r)] + list(re) + [numpy.interp(sr_top, s, r)]
    else:
        snew = [sr_base, sr_top]
        rnew = [numpy.interp(sr_base, s, r), numpy.interp(sr_top, s, r)]

    S = trapz(rnew, snew) * Lshape * Lwshape

    return S


def form_factor(leaf):
    _, _, s, r = leaf
    return blade_elt_area(s, r, 1, 1, 0, 1)


def arrange_leaf(leaf, stem_diameter=0, inclination=1, relative=True):
    """Arrange a leaf to be placed along a stem with a given inclination.

    Args:
        leaf: a x, y, s, r tuple describing leaf shape
        stem_diameter: the diameter of the sem at the leaf insertion point
        inclination: if relative=False, the leaf basal inclination (deg). A
        multiplier to leaf basal inclination angle otherwise
        relative: (bool) controls the meaning of inclination parameter

    Returns:
        a modified x, y, s, r tuple

    """

    x, y, s, r = map(numpy.array, leaf)
    if relative and inclination == 1:
        x1, y1 = x, y
    else:
        basal_inclination = pgl.angle((x[1] - x[0], y[1] - y[0]), (0, 1))

        if relative:
            angle = inclination * basal_inclination
            angle = min(pi, angle)
        else:
            angle = radians(inclination)

        rotation_angle = basal_inclination - angle

        # rotation of the midrib
        cos_a = cos(rotation_angle)
        sin_a = sin(rotation_angle)

        x1 = x[0] + cos_a * x - sin_a * y
        y1 = y[0] + sin_a * x + cos_a * y
    leaf = x1 + stem_diameter / 2., y1, s, r

    return leaf

def leaf_mesh(leaf, L_shape, Lw_shape, length, s_base=0, s_top=1, flipx=False,
              twist=0, volume=0, stem_diameter=0, inclination=1,
              relative=False):
    """Compute mesh for a leaf element along a scaled leaf shape

    Args:
        leaf: a x,y,s,r tuple describing leaf shape
        L_shape: length of the shape
        Lw_shape: width of the shape
        length: the total visible length to be meshed
        s_base: normalised position (on the shape) of the start of the element
        s_top: normalised position (on the shape) of the end of the element
        flipx: should leaf shape be flipped ?
        twist:
        volume: float value of the tickness of the leaf.
              Default is 0. Else it indicates the relative depth of the leaf
              along the midrib.
        stem_diameter: the diameter of the sem at the leaf insertion point
        inclination: if relative=False, the leaf basal inclination (deg). A
        multiplier to leaf basal inclination angle otherwise
        relative: (bool) controls the meaning of inclination parameter

    Returns:
        a PlanGl mesh representing the element
    """
    shape = arrange_leaf(leaf, stem_diameter=float(stem_diameter) / L_shape,
                         inclination=inclination, relative=relative)
    # flip to position leaves along tiller emitted
    if flipx:
        # to position leaves along tiller emitted
        shape = (-shape[0],) + shape[1:]

    mesh = fitting.mesh4(shape,
                         L_shape,
                         length,
                         s_base,
                         s_top,
                         Lw_shape,
                         twist=twist,
                         volume=volume)

    if mesh:
        pts, ind = mesh
        if len(ind) < 1:
            mesh = None
        else:
            mesh = fitting.plantgl_shape(pts, ind)
    else:
        if length > 0:
            print('ERROR No mesh', s_base, s_top, length)
            pass
        mesh = None

    return mesh