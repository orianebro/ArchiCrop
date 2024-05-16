"""Generate a geometric-based MTG representation of a cereal plant"""
from __future__ import annotations

import numpy as np
import pandas as pd
from openalea.mtg import MTG, fat_mtg
from scipy.interpolate import interp1d

from .geometry import mtg_interpreter
from .plant_design import get_form_factor


def curvilinear_abscisse(x, y, z=None):
    """Curvilinear abcissa along a polyline"""
    s = np.zeros(len(x))
    if z is None:
        s[1:] = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
    else:
        s[1:] = np.sqrt(
            np.diff(x) ** 2 + np.diff(y) ** 2 + np.diff(z) ** 2)
    return s.cumsum()


def majors_axes_regression(x, y):
    """ Performs a major axis regression
    :return: a, b, c : (float) coef of the regression line ax + by + c = 0
    """

    x = np.array(x)
    y = np.array(y)
    xm = np.mean(x)
    ym = np.mean(y)
    s_xy = ((x - xm) * (y - ym)).sum()
    s_xx = np.power(x - xm, 2).sum()
    s_yy = np.power(y - ym, 2).sum()

    if s_xx == 0:
        a = 1
        b = 0
        c = -xm
    else:
        b = -1
        a = np.sqrt(s_yy / s_xx) if s_xy > 0 else -np.sqrt(s_yy / s_xx)
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

    x, y, s, r = list(map(np.array, leaf))
    cx = curvilinear_abscisse(x, y)
    cx_m = max(cx)
    xo, yo, zo = origin
    azimuth = np.radians(azimuth)
    leaf_x = xo + x / cx_m * length * np.cos(azimuth)
    leaf_y = yo + x / cx_m * length * np.sin(azimuth)
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

    x, y, z, r = list(map(lambda x : np.array(x).astype(float), list(zip(*polyline))))
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
    azimuth = np.degrees(np.arctan2(yp[-1] - yp[0], xp[-1] - xp[0]))
    return (x_leaf / sxp.max(), y_leaf / sxp.max(), s, r), length, radius_max, azimuth, origin


def as_json(plant):
    """convert a plant dimension + shape table in a cereals json input
    (reverse of 'as_plant')"""
    internode = plant.L_internode.values
    diameter = plant.W_internode.values
    stem = [0, *internode.cumsum().tolist()]
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
    df = pd.DataFrame(
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
    stem = [0, *df.hins.tolist()]
    df['internode'] = np.diff(stem)
    df['ntop'] = df.index.max() - df.index + 1
    # re-index leaves with ntop
    leaves = {df.ntop[rank] : leaves[rank] for rank in df.index}
    blades = pd.DataFrame({'L_blade': df.l_leaf,
                               'S_blade': df.area,
                               'W_blade': df.w_leaf,
                               'ntop': df.ntop,
                               'plant': 1,
                               'leaf_azimuth': df.azimuth,
                               'form_factor' : df.ff})

    stem = pd.DataFrame({'L_internode': df.internode,
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
        relative_azimuth[1:] = np.diff(relative_azimuth)
    else:
        dim = plant
        leaves = {row['ntop']: row['leaf_shape'] for index, row in dim.iterrows()}

    # print(dim)

    g = MTG()
    vid_plant = g.add_component(g.root, label='Plant', edge_type='/') 
    vid_axis = g.add_component(vid_plant, label='MainAxis', edge_type='/') 


    first_internode = True

    for _i, row in dim.iterrows():
        internode = {'label': 'Stem',
                     'mature_length': row['L_internode'],
                     'length': row['L_internode'],
                     'visible_length': row['L_internode'],
                     'is_green': True,
                     'diameter': row['W_internode'],
                     'azimuth': row['leaf_azimuth'],
                     'grow': False}

        if first_internode :
            # vid_metamer = g.add_component(vid_axis)
            vid_internode = g.add_component(vid_axis, **internode)
            # g.node(vid_metamer).label='Metamer'
            # g.node(vid_metamer).edge_type='/'
            first_internode=False
        else:
            vid_internode = g.add_child(vid_internode, edge_type='<', **internode)
            # vid_internode, vid_metamer = g.add_child_and_complex(vid_internode, edge_type='<', **internode)
            # g.node(vid_metamer).label='Metamer'
            # g.node(vid_metamer).edge_type='<'

        leaf = {'label': 'Leaf',
                'shape': leaves[row['ntop']],
                'shape_mature_length': row['L_blade'],
                'length': row['L_blade'],
                # 'visible_length': row['L_blade'],
                'visible_length': row['L_blade'],
                'leaf_area' : row['S_blade'],
                'form_factor' : row['form_factor'],
                'is_green': True,
                'srb': 0,
                'srt': 1,
                'lrolled': 0,
                'd_rolled': 0,
                'shape_max_width': row['W_blade'],
                'stem_diameter': row['W_internode'],
                'grow': False}

        vid_leaf = g.add_child(vid_internode, edge_type='+', **leaf)  # noqa: F841

    g = fat_mtg(g)

    # Compute geometry
    g = mtg_interpreter(g, classic=classic)

    return g  # noqa: RET504


## from adel, for inspiration

'''
def mtg_factory(parameters, metamer_factory=adel_metamer, leaf_sectors=1,
                leaves=None, stand=None, axis_dynamics=None,
                add_elongation=False, topology=['plant', 'axe_id', 'numphy'],
                split=False, aborting_tiller_reduction=1.0, leaf_db=None):
    """ Construct a MTG from a dictionary of parameters.

    The dictionary contains the parameters of all metamers in the stand (topology + properties).
    metamer_factory is a function that build metamer properties and metamer elements from parameters dict.
    leaf_sectors is an integer giving the number of LeafElements per Leaf blade
    leaves is a {species:adel.geometric_elements.Leaves} dict
    stand is a list of tuple (xy_position_tuple, azimuth) of plant positions
    axis_dynamics is a 3 levels dict describing axis dynamic. 1st key level is plant number, 2nd key level is axis number, and third ky level are labels of values (n, tip, ssi, disp)
    topology is the list of key names used in parameters dict for plant number, axe number and metamer number
    aborting_tiller_reduction is a scaling factor applied to reduce all dimensions of organs of tillers that will abort

    Axe number 0 is compulsory
  
    """

    def get_component(components, index):
        component = components[index]
        elements = component['elements']
        properties = dict(component)
        del properties['elements']
        return properties, elements

    if leaves is None:
        dynamic_leaf_db = {0: False}
        leaves = {0: None}
    else:
        dynamic_leaf_db = {k: leaves[k].dynamic for k in leaves}

    g = MTG()

    # buffers
    # for detection of newplant/newaxe
    prev_plant = 0
    prev_axe = -1
    # current vid
    vid_plant = -1
    vid_axe = -1
    vid_metamer = -1
    vid_node = -1
    vid_elt = -1
    # vid of top of stem nodes and elements
    vid_topstem_node = -1
    vid_topstem_element = -1
    # vid of plant main stem (axe0)
    vid_main_stem = -1
    # buffer for the vid of main stem anchor points for the first metamer, node and element of tillers
    metamers = []
    nodes = []
    elts = []

    dp = parameters
    nrow = len(dp['plant'])

    for i in range(nrow):
        plant, num_metamer = [int(convert(dp.get(x)[i], undef=None)) for x in
                              [topology[e] for e in [0, 2]]]
        axe = dp.get(topology[1])[i]
        mspos = int(convert(dp.get('ms_insertion')[i], undef=None))
        args = properties_from_dict(dp, i, exclude=topology)
        # Add plant if new
        if plant != prev_plant:
            label = 'plant' + str(plant)
            position = (0, 0, 0)
            azimuth = 0
            species = 0
            if 'species' in args:
                species = args['species']
            if stand and len(stand) >= plant:
                position, azimuth = stand[plant - 1]
            vid_plant = g.add_component(g.root, label=label, edge_type='/',
                                        position=position, azimuth=azimuth,
                                        refplant_id=args.get('refplant_id'),
                                        species=species)
            # reset buffers
            prev_axe = -1
            vid_axe = -1
            vid_metamer = -1
            vid_node = -1
            vid_elt = -1
            vid_topstem_node = -1
            vid_topstem_element = -1
            vid_main_stem = -1
            metamers = []
            nodes = []
            elts = []

        # Add axis
        if axe != prev_axe:
            label = ''.join(axe.split('.'))
            timetable = None
            if axis_dynamics:
                timetable = axis_dynamics[str(plant)][str(axe)]
            if axe == 'MS':
                vid_axe = g.add_component(vid_plant, edge_type='/', label=label,
                                          timetable=timetable,
                                          HS_final=args.get('HS_final'),
                                          nff=args.get('nff'),
                                          hasEar=args.get('hasEar'),
                                          azimuth=args.get('az_insertion'))
                vid_main_stem = vid_axe
            else:
                vid_axe = g.add_child(vid_main_stem, edge_type='+', label=label,
                                      timetable=timetable,
                                      HS_final=args.get('HS_final'),
                                      nff=args.get('nff'),
                                      hasEar=args.get('hasEar'),
                                      azimuth=args.get('az_insertion'))

        # Add metamer
        assert num_metamer > 0
        # args are added to metamers only if metamer_factory is none, otherwise compute metamer components
        components = []
        if metamer_factory:
            xysr_key = None
            if leaves[species] is not None and 'LcType' in args and 'LcIndex' in args:
                lctype = int(args['LcType'])
                lcindex = int(args['LcIndex'])
                if lctype != -999 and lcindex != -999:
                    age = None
                    if dynamic_leaf_db[species]:
                        age = float(args[
                                        'rph']) - 0.3  # age_db = HS - rank + 1 = ph - 1.3 - rank +1 = rph - .3
                        if age != 'NA':
                            age = max(0, int(float(age)))
                    xysr_key = leaves[species].get_leaf_key(lctype, lcindex,
                                                            age)

            elongation = None
            if add_elongation:
                startleaf = -.4
                endleaf = 1.6
                stemleaf = 1.2
                startE = endleaf
                endE = startE + (endleaf - startleaf) / stemleaf
                endBlade = endleaf
                if args['Gl'] > 0:
                    endBlade = args['Ll'] / args['Gl'] * (endleaf - startleaf)
                elongation = {'startleaf': startleaf, 'endBlade': endBlade,
                              'endleaf': endleaf, 'endE': endE}
            if not 'ntop' in args:
                args.update({'ntop': None})
            if not 'Gd' in args:
                args.update({'Gd': 0.19})
            args.update({'split': split})
            
            hs_f = args.get('HS_final')
            if hs_f != 'NA':
                if float(hs_f) < args.get('nff'):
                    for what in (
                    'Ll', 'Lv', 'Lr', 'Lsen', 'L_shape', 'Lw_shape', 'Gl', 'Gv',
                    'Gsen', 'Gd', 'El', 'Ev', 'Esen', 'Ed'):
                        args.update(
                            {what: args.get(what) * aborting_tiller_reduction})
            components = metamer_factory(Lsect=leaf_sectors, shape_key=xysr_key,
                                         elongation=elongation,
                                         leaves=leaves[species], **args)
            args = {'L_shape': args.get('L_shape')}
        #
        label = 'metamer' + str(num_metamer)
        new_metamer = g.add_component(vid_axe, edge_type='/', label=label,
                                      **args)
        if axe == 'MS' and num_metamer == 1:
            vid_metamer = new_metamer
        elif num_metamer == 1:
            # add the edge with the bearing metamer on main stem
            vid_metamer = metamers[mspos - 1]
            vid_metamer = g.add_child(vid_metamer, child=new_metamer,
                                      edge_type='+')
        else:
            vid_metamer = g.add_child(vid_metamer, child=new_metamer,
                                      edge_type='<')

        # add metamer components, if any           
        if len(components) > 0:
            # deals with first component (internode) and first element 
            node, elements = get_component(components, 0)
            element = elements[0]
            new_node = g.add_component(vid_metamer, edge_type='/', **node)
            new_elt = g.add_component(new_node, edge_type='/', **element)
            if axe == 'MS' and num_metamer == 1:  # root of main stem
                vid_node = new_node
                vid_elt = new_elt
            elif num_metamer == 1:  # root of tiller
                vid_node = nodes[mspos - 1]
                vid_node = g.add_child(vid_node, child=new_node, edge_type='+')
                vid_elt = elts[mspos - 1]
                vid_elt = g.add_child(vid_elt, child=new_elt, edge_type='+')
            else:
                vid_node = g.add_child(vid_topstem_node, child=new_node,
                                       edge_type='<')
                vid_elt = g.add_child(vid_topstem_element, child=new_elt,
                                      edge_type='<')
            # add other elements of first component (the internode)
            for i in range(1, len(elements)):
                element = elements[i]
                vid_elt = g.add_child(vid_elt, edge_type='<', **element)
            vid_topstem_node = vid_node
            vid_topstem_element = vid_elt  # last element of internode

            # add other components   
            for i in range(1, len(components)):
                node, elements = get_component(components, i)
                if node['label'] == 'sheath':
                    edge_type = '+'
                else:
                    edge_type = '<'
                vid_node = g.add_child(vid_node, edge_type=edge_type, **node)
                element = elements[0]
                new_elt = g.add_component(vid_node, edge_type='/', **element)
                vid_elt = g.add_child(vid_elt, child=new_elt,
                                      edge_type=edge_type)
                for j in range(1, len(elements)):
                    element = elements[j]
                    vid_elt = g.add_child(vid_elt, edge_type='<', **element)

                    # update buffers
        if axe == 'MS':
            metamers.append(vid_metamer)
            if len(components) > 0:
                nodes.append(vid_topstem_node)
                elts.append(vid_topstem_element)
        prev_plant = plant
        prev_axe = axe

    return fat_mtg(g)
'''