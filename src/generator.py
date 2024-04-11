"""Generate a geometric-based MTG representation of a cereal plant"""
import numpy
import pandas

from scipy.interpolate import interp1d
from openalea.mtg import MTG, fat_mtg
from geometry import mtg_interpreter
from plant_design import get_form_factor


def curvilinear_abscisse(x, y, z=None):
    """Curvilinear abcissa along a polyline"""
    s = numpy.zeros(len(x))
    if z is None:
        s[1:] = numpy.sqrt(numpy.diff(x) ** 2 + numpy.diff(y) ** 2)
    else:
        s[1:] = numpy.sqrt(
            numpy.diff(x) ** 2 + numpy.diff(y) ** 2 + numpy.diff(z) ** 2)
    return s.cumsum()


def majors_axes_regression(x, y):
    """ Performs a major axis regression
    :return: a, b, c : (float) coef of the regression line ax + by + c = 0
    """

    x = numpy.array(x)
    y = numpy.array(y)
    xm = numpy.mean(x)
    ym = numpy.mean(y)
    s_xy = ((x - xm) * (y - ym)).sum()
    s_xx = numpy.power(x - xm, 2).sum()
    s_yy = numpy.power(y - ym, 2).sum()

    if s_xx == 0:
        a = 1
        b = 0
        c = -xm
    else:
        b = -1
        if s_xy > 0:
            a = numpy.sqrt(s_yy / s_xx)
        else:
            a = -numpy.sqrt(s_yy / s_xx)
        c = ym - a * xm
    return a, b, c


def line_projection(a, b, c, xo, yo):
    """coordinate of the projection of xo, yo on ax + by + c = 0 line"""
    x = (b * (b * xo - a * yo) - a * c) / (a**2 + b**2)
    y = (a * (-b * xo + a* yo) - b * c) / (a**2 + b**2)
    return x, y


def as_polyline(leaf, length=1, radius_max=1, origin=(0, 0, 0), azimuth=0):
    """Transform x, y, s, r leaf tuple into a [(x,y,z,r), ...] polyline
    (reverse of 'as_leaf')

    Args:
        leaf: a (x, y, s, r) tuple giving x, y coordinates of leaf
        midrib in the frame defined by vertical leaf plane and leaf insertion
        point as origin. s, r are normalised curvilinear abscissa and radius
        along midrib
        length: leangth of the leaf
        radius_max:  maximal width of the leaf
        origin: the coordinate of leaf insertion point
        azimuth: the angle (deg, from X+ positive counter-clockwise)

    Returns:
        a (x,y, z,r) list of tuple sampling leaf polyline

    """

    x, y, s, r = list(map(numpy.array, leaf))
    cx = curvilinear_abscisse(x, y)
    cx_m = max(cx)
    xo, yo, zo = origin
    azimuth = numpy.radians(azimuth)
    leaf_x = xo + x / cx_m * length * numpy.cos(azimuth)
    leaf_y = yo + x / cx_m * length * numpy.sin(azimuth)
    leaf_z = zo + y / cx_m * length
    leaf_r = r * radius_max
    return list(zip(*list(map(tuple, (leaf_x,
                            leaf_y,
                            leaf_z,
                            leaf_r)))))


def as_leaf(polyline):
    """Compute leaf x,y,s,r tuple, length, radius max and azimuth from a
    polyline (reverse of 'as_polyline')

    Args:
        polyline: a [(x,y,z,r), ...] list of tuple sampling leaf midrib

    Returns:
        (x, y, s, r), length, radius_max , azimuth, origin of the leaf
    """

    x, y, z, r = list(map(lambda x : numpy.array(x).astype(float), list(zip(*polyline))))
    xo, yo, zo = x[0], y[0], z[0]
    sx = curvilinear_abscisse(x, y, z)
    a, b, c = majors_axes_regression(x, y)
    #
    length = sx.max()
    radius_max = r.max()
    s = sx / length
    r /= radius_max
    origin = (xo, yo, zo)
    y_leaf = z - zo
    xp, yp = list(zip(*[line_projection(a, b, c, x[0], x[1]) for x in zip(x,y)]))
    x_leaf = curvilinear_abscisse(xp, yp)
    sxp = curvilinear_abscisse(x_leaf, y_leaf)
    azimuth = numpy.degrees(numpy.arctan2(yp[-1] - yp[0], xp[-1] - xp[0]))
    return (x_leaf / sxp.max(), y_leaf / sxp.max(), s, r), length, radius_max, azimuth, origin


def as_json(plant):
    """convert a plant dimension + shape table in a cereals json input
    (reverse of 'as_plant')"""
    internode = plant.L_internode.values
    diameter = plant.W_internode.values
    stem = [0] + internode.cumsum().tolist()
    stem_diameter = [diameter[0]] + diameter.tolist()
    polylines = [as_polyline(leaf, length, width, (0,0,h), azim) for leaf, length, width, h, azim
                 in
                 zip(plant.leaf_shape, plant.L_blade, plant.W_blade,
                     plant.h_ins, plant.leaf_azimuth)]

    return {
        'leaf_polylines': polylines,
        'leaf_order': plant.leaf_rank.values.tolist(),
        'stem': [(0, 0, z, r) for z, r in zip(stem, stem_diameter)]}


def as_plant(json):
    """restore plant dimension + leaves representation of a plant encoded in
    json (reverse of as_json)"""

    ranks = json['leaf_order']
    leaves, l_leaf, w_leaf, azimuth, origin = list(zip(*list(map(as_leaf, json['leaf_polylines']))))
    # use rank as index
    df = pandas.DataFrame(
        {'rank': ranks,
         'l_leaf': l_leaf,
         'w_leaf': w_leaf,
         'azimuth': azimuth,
         'hins': [ori[2] for ori in origin]}).set_index('rank')
    df.sort_index(inplace=True)
    leaves = {rank:leaves[i] for i, rank in enumerate(ranks)}
    x_stem, y_stem, z_stem, r_stem = list(zip(*json['stem']))
    df['diam'] = interp1d(
        z_stem, r_stem, bounds_error=False, fill_value=r_stem[-1])(df.hins)
    df['ff'] = [get_form_factor(leaves[r]) for r in df.index]
    df['area'] = df.l_leaf * df.w_leaf * df.ff
    stem = [0] + df.hins.tolist()
    df['internode'] = numpy.diff(stem)
    df['ntop'] = df.index.max() - df.index + 1
    # re-index leaves with ntop
    leaves = {df.ntop[rank] : leaves[rank] for rank in df.index}
    blades = pandas.DataFrame({'L_blade': df.l_leaf,
                               'S_blade': df.area,
                               'W_blade': df.w_leaf,
                               'ntop': df.ntop,
                               'plant': 1,
                               'leaf_azimuth': df.azimuth})

    stem = pandas.DataFrame({'L_internode': df.internode,
                             'L_sheath': 0,
                             'W_internode': df.diam,
                             'W_sheath': df.diam,
                             'h_ins': df.hins,
                             'ntop': df.ntop,
                             'plant': 1})
    return blades, stem, leaves


def cereals(json=None,
            classic=False,
            unit='cm',
            seed=None,
            plant=None):
    """
    Generate a 'geometric-based' MTG representation of cereals

    Args:
        json: a dict of parameters with:
            leaf_order : a list of int defining leaf rank
            leaf_polylines : [[(x, y, z, r), ..., ], ...] list of list tuple
            sampling  leaf midribs. r is leaf width at position x, y, z along
            leaf midrib
            stem : [(x, y, z, r), ..., ] list of tuples sampling stem from its
            base to its end. r is the stem diameter at x,y,z position along the
             stem
        classic: (bool) should stem cylinders be classical pgl cylinders ?
        unit: (string) desired length unit for the output mtg
        seed: (int) a seed for the random number generator
    """
    if plant is None:
        blade_dimensions, stem_dimensions, leaves = as_plant(json)

        dim = blade_dimensions.merge(stem_dimensions)
        dim = dim.sort_values('ntop', ascending=False)
        relative_azimuth = dim.leaf_azimuth.copy()
        relative_azimuth[1:] = numpy.diff(relative_azimuth)
    else:
        dim = plant
        leaves = {row['ntop']: row['leaf_shape'] for index, row in dim.iterrows()}


    g = MTG()
    vid_node = g.add_component(g.root, label='plant', edge_type='/')
    for i, row in dim.iterrows():
        internode = {'label': 'StemElement',
                     'length': row['L_internode'],
                     'is_green': True,
                     'diameter_base': row['W_internode'],
                     'diameter_top': row['W_internode'],
                     'azimuth': row['leaf_azimuth']}

        vid_node = g.add_child(vid_node, edge_type='<', **internode)

        blade = {'label': 'LeafElement',
                 'shape': leaves[row['ntop']],
                 'shape_mature_length': row['L_blade'],
                 'length': row['L_blade'],
                 'visible_length': row['L_blade'],
                 'is_green': True,
                 'srb': 0,
                 'srt': 1,
                 'lrolled': 0,
                 'd_rolled': 0,
                 'shape_max_width': row['W_blade'],
                 'stem_diameter': row['W_internode']}

        vid_blade = g.add_child(vid_node, edge_type='+', **blade)

    g = fat_mtg(g)

    # Compute geometry
    g = mtg_interpreter(g, classic=classic)

    return g
