"""Generate a geometric-based MTG representation of a cereal plant"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from openalea.mtg import MTG, fat_mtg

from .cereal_leaf import (
    blade_dimension,
    get_form_factor,
    leaf_as_dict,
    parametric_leaf,
    shape_to_surface,
)
from .cereal_plant import bell_shaped_dist, geometric_dist, leaf_azimuth
from .fitting import curvilinear_abscisse
from .geometry import mtg_interpreter
from .internode import stem_as_dict, stem_dimension


def stem_element_properties(nb_phy, nb_short_phy, short_phy_height, height, stem_q, diam_top, diam_base, ntop):
    '''
    Return dict of stem element properties
    '''
    # short_phy_height = 2
    pseudostem_height = nb_short_phy * short_phy_height

    pseudostem = np.array([short_phy_height*i for i in range(1, nb_short_phy+1)])
    stem = np.array(geometric_dist(height, nb_phy-nb_short_phy, q=stem_q, u0=pseudostem_height))
    insertion_heights = np.concatenate((pseudostem, stem), axis=0)

    stem_diameters = np.linspace(diam_base, diam_top, nb_phy).tolist()

    return stem_dimension(
        h_ins=insertion_heights, d_internode=stem_diameters, ntop=ntop
    )


def leaf_properties(nb_phy, leaf_area, rmax, skew, insertion_angle, scurv, curvature, 
                    klig, swmax, f1, f2, ntop, wl, phyllotactic_angle, phyllotactic_deviation, plant_orientation, spiral):
    '''
    Return dict of blade properties
    '''
    # leaf areas
    leaf_areas = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))

    # leaf shapes
    a_leaf = parametric_leaf(
        nb_segment=10,
        insertion_angle=insertion_angle,
        scurv=scurv,
        curvature=curvature,
        klig=klig, swmax=swmax, f1=f1, f2=f2
    )

    blades = blade_dimension(area=leaf_areas, ntop=ntop, wl=wl)

    # phyllotaxy
    leaf_azimuths = leaf_azimuth(
        size=nb_phy,
        phyllotactic_angle=phyllotactic_angle,
        phyllotactic_deviation=phyllotactic_deviation,
        plant_orientation=plant_orientation,
        spiral=spiral,
    )
    leaf_azimuths[1:] = np.diff(leaf_azimuths)

    blades["leaf_azimuth"] = leaf_azimuths
    blades["leaf_shape"] = [a_leaf] * nb_phy 
    blades["tck"] = [shape_to_surface(shape, wl) for shape in blades["leaf_shape"]]

    return blades


def add_development(g, vid, tt, dtt, rate):
    """
    Add dynamic properties on the mtg to simulate development

    :param g: MTG, MTG of a plant
    :param vid: int, vertex id of the plant component
    :param tt: float, current thermal time (in °C.day)
    :param dtt: float, thermal time increment (in °C.day)
    :param rate: float, rate of development (in °C.day)

    :return: tt
    """
    g.node(vid).start_tt = tt
    g.node(vid).end_tt = tt + dtt
    tt += rate
    return tt


def add_leaf_senescence(g, vid_leaf, leaf_lifespan, end_juv):
    if isinstance(leaf_lifespan, list): 
        if g.node(vid_leaf).start_tt < end_juv:
            g.node(vid_leaf).senescence = g.node(vid_leaf).start_tt + leaf_lifespan[0] 
        else:
            g.node(vid_leaf).senescence = g.node(vid_leaf).start_tt + leaf_lifespan[1]
    else:
        g.node(vid_leaf).senescence = g.node(vid_leaf).start_tt + leaf_lifespan



def add_tiller(g, vid, start_time, phyllochron, plastochron, 
               stem_duration, leaf_duration, leaf_lifespan, end_juv, 
               tiller_delay, reduction_factor,  # noqa: ARG001
               height, leaf_area, nb_short_phy, short_phy_height, wl, diam_base, diam_top,
               insertion_angle, scurv, curvature,
               klig, swmax, f1, f2,
               stem_q, rmax, skew, 
               phyllotactic_angle, phyllotactic_deviation,
               tiller_angle, gravitropism_coefficient,
               plant_orientation=45,
               spiral=True):
    """ Add a tiller to the plant at the given time
    Args:
        g: MTG
        vid: vertex id of the plant
        start_time: time of tiller initiation (in °C.day)
        phyllochron: phyllochron of the plant (in °C.day/stem element)
        tiller_delay: delay of the tiller (in °C.day)
        plastochron: plastochron of the plant (in °C.day/leaf)
        stem_duration: duration of the stem element development (in phyllochron)
        leaf_duration: duration of the leaf blade development (in phyllochron)
        leaf_lifespan: lifespan of the leaf (in °C.day)
        end_juv: end of the juvenile phase (in °C.day)
        reduction_factor: factor of reduction for tiller properties
        height: height of the plant (in cm)
        leaf_area: maximal potential leaf area of the plant (in cm²)
        nb_short_phy: number of short phytomers
        short_phy_height: height of the short phytomers (in cm)
        wl: leaf width-to-length ratio
        diam_base: diameter of the base of the main stem (in cm)
        diam_top: diameter of the top of the main stem (in cm)
        insertion_angle: insertion angle of the leaf (i.e. between the stem and the tangent line at the base of the leaf) (in °)
        scurv: curvilinear abscissa of inflexion point for leaf curvature (in [0,1])
        curvature: curvature angle (i.e. angle between insertion angle and the tangent line at the tip of the leaf) (in °)
        klig: coefficient for leaf shape
        swmax: relative leaf length where leaf reaches its maximum width (in [0,1])
        f1: coefficient for leaf shape
        f2: coefficient for leaf shape
        stem_q: common ratio of the geometric series defining stem element lengths for a given height (cf partition of unit)
        rmax: relative leaf rank corresponding to the position of the longest leaf, from the base of the stem, according to a bell-shaped distribution of leaf lengths along the stem (in [0,1])
        skew: parameter describing the asymmetry of the bell-shaped distribution of leaf lengths along the stem
        phyllotactic_angle: angle between the midribs of two consecutive leaves around the stem (in °)
        phyllotactic_deviation: half-amplitude of deviation around phyllotactic angle (in °)
        tiller_angle: angle of the tiller wrt to tiller of prior order (in °)
        gravitropism_coefficient: coefficient of gravitropism for the tiller
        plant_orientation: orientation of the plant (in °)
        spiral: whether the phyllotaxy is spiral or not
        
    Returns:
        List[Tuple[int, float]]: List of tuples containing the vertex id of the tiller and its start time
    """

    # Add a new component to the MTG

    tillers = []
    # scale_id = g.scale(vid)
    axis_id = g.complex(vid)
    rank = g.Rank(vid) + 1  # Number of edges from the root of the axis
    n = len(g.Axis(vid))
    len_tiller  = n - rank  # we remove the parent that do not belong to the tiller 
    nb_phy = len_tiller

    ranks = range(1, nb_phy + 1)
    ntop = max(ranks) - np.array(ranks) + 1

    # see how properties change with reduction factor and order p.r^o
    stem_prop = stem_element_properties(nb_phy, nb_short_phy, short_phy_height, height, stem_q, diam_top, diam_base, ntop)
    leaf_prop = leaf_properties(nb_phy, leaf_area, rmax, skew, insertion_angle, scurv, curvature, 
                    klig, swmax, f1, f2, ntop, wl, phyllotactic_angle, phyllotactic_deviation, plant_orientation, spiral)

    tt_stem = start_time
    tt_leaf = start_time
    dtt_stem = phyllochron * stem_duration
    dtt_leaf = plastochron * leaf_duration

    tid = g.add_child(parent=axis_id, edge_type='+', label='Axis')

    first = True

    for rank in ranks:

        stem = stem_as_dict(stem_prop=stem_prop, leaf_prop=leaf_prop, rank=rank)
        if first:
            vid_stem, tid2 = g.add_child_and_complex(parent=vid, complex=tid, edge_type='+', **stem)
            g.node(vid_stem).tiller_angle = tiller_angle
            g.node(vid_stem).gravitropism_coefficient = gravitropism_coefficient
            first = False
        else:
            vid_stem = g.add_child(parent=vid_stem, edge_type='<', **stem)
        tt_stem = add_development(g=g, vid=vid_stem, tt=tt_stem, dtt=dtt_stem, rate=phyllochron)

        leaf = leaf_as_dict(stem_prop=stem_prop, leaf_prop=leaf_prop, rank=rank, wl=wl)
        vid_leaf = g.add_child(vid_stem, edge_type="+", **leaf)
        tt_leaf = add_development(g=g, vid=vid_leaf, tt=tt_leaf, dtt=dtt_leaf, rate=plastochron)
        add_leaf_senescence(g=g, vid_leaf=vid_leaf, leaf_lifespan=leaf_lifespan, end_juv=end_juv)
        
        tillers.append((vid_stem, tt_stem + tiller_delay))

    return tillers


def cereals(nb_phy, phyllochron, plastochron, stem_duration, leaf_duration, 
            leaf_lifespan, end_juv, nb_tillers, tiller_delay, reduction_factor,
            height,
            leaf_area,
            nb_short_phy,
            short_phy_height,
            wl,
            diam_base,
            diam_top,
            insertion_angle,
            scurv,
            curvature,
            # alpha,
            klig, swmax, f1, f2,
            stem_q,
            rmax,
            skew,
            phyllotactic_angle,
            phyllotactic_deviation,
            tiller_angle,
            gravitropism_coefficient,
            plant_orientation=45,
            spiral=True,
            classic=False):
            # json=None, classic=False, seed=None, plant=None):
    """
    Generate a 'geometric-based' MTG representation of cereals

    Args:
    :param nb_phy: int, number of phytomers
    :param height: float, height of the main stem (in cm)
    :param leaf_area: float, maximal potential leaf area of plant (in cm²)
    :param wl: float, leaf width-to-length ratio
    :param diam_base: float, diameter of the base of the main stem (in cm)
    :param diam_top: float, diameter of the top of the main stem (in cm)
    :param insertion_angle: float, insertion angle of the leaf (i.e. between the stem and the tangent line at the base of the leaf) (in °)
    :param scurv: float, curvilinear abscissa of inflexion point for leaf curvature (in [0,1])
    :param curvature: float, curvature angle (i.e. angle between insertion angle and the tangent line at the tip of the leaf) (in °)
    :param stem_q: float, common ratio of the geometric series defining internode lengths for a given height (cf partition of unit)
    :param rmax: float, proportion of the total number of leaves corresponding to the position of the longest leaf, from the base of the stem, according to a bell-shaped distribution of leaf lengths along the stem (in [0,1])
    :param skem: float, parameter describing the asymmetry of the bell-shaped distribution of leaf lengths along the stem
    :param phyllotactic_angle: float, angle between the midribs of two consecutive leaves around the stem (in °)
    :param phyllotactic_deviation: float, half-amplitude of deviation around phyllotactic angle (in °)
    :param phyllochron: float, phyllochron, i.e. internode appearance rate (in °C.day/internode)
    :param plastochron: float, plastochron, i.e. leaf appearance rate (in °C.day/leaf)
    :param leaf_duration: float, phyllochronic time for a leaf to develop from tip appearance to collar appearance (/phyllochron)
    :param stem_duration: float, phyllochronic time for a stem to develop from base to top (/phyllochron)
    :param leaf_lifespan: float, lifespan of a leaf (in °C.day)
    :param end_juv: float, end of juvenile phase (in °C.day)
    :param nb_tillers: int, number of tillers 
    :param tiller_delay: float, delay between the appearance of a phytomer and the appearance of a tiller from its lateral meristem (/phyllochron)
    :param reduction_factor: float, factor of reduction for tiller properties regarding the main stem properties
    :param plant_orientation: float, orientation of the plant (in °)
    :param spiral: bool, whether the phyllotaxy is spiral or not
    :param classic: bool, whether to use the classic MTG interpreter or not

    Returns:
    :return: MTG, a geometric-based MTG representation of cereals
    """

    # Main Axis

    ranks = range(1, nb_phy + 1)
    ntop = max(ranks) - np.array(ranks) + 1
    # Dicts instead of a dataframe
    stem_prop = stem_element_properties(nb_phy=nb_phy, nb_short_phy=nb_short_phy, short_phy_height=short_phy_height, height=height, 
                                        stem_q=stem_q, diam_top=diam_top, diam_base=diam_base, ntop=ntop)
    leaf_prop = leaf_properties(nb_phy=nb_phy, leaf_area=leaf_area, rmax=rmax, skew=skew, insertion_angle=insertion_angle, 
                                scurv=scurv, curvature=curvature, klig=klig, swmax=swmax, f1=f1, f2=f2, ntop=ntop, wl=wl, 
                                phyllotactic_angle=phyllotactic_angle, phyllotactic_deviation=phyllotactic_deviation, 
                                plant_orientation=plant_orientation, spiral=spiral)
    
    tt_stem = 0
    tt_leaf = 0
    dtt_stem = phyllochron * stem_duration
    dtt_leaf = plastochron * leaf_duration

    tiller_points = []

    g = MTG()
    # Add a root vertex for the plant
    vid_plant = g.add_component(g.root, label="Plant", edge_type="/") 
    # Add a plant vertex for the main axis
    vid_axis = g.add_component(vid_plant, label="MainAxis", edge_type="/")

    first = True

    # iterate over the number of phytomers of the main stem
    for rank in ranks:

        stem = stem_as_dict(stem_prop=stem_prop, leaf_prop=leaf_prop, rank=rank)
        if first:
            vid_stem = g.add_component(vid_axis, **stem)
            g.node(vid_stem).tiller_angle = 0
            g.node(vid_stem).gravitropism_coefficient = 0
            first = False
        else:
            vid_stem = g.add_child(vid_stem, edge_type="<", **stem)
        tt_stem = add_development(g=g, vid=vid_stem, tt=tt_stem, dtt=dtt_stem, rate=phyllochron)
        tiller_points.append((vid_stem, tt_stem + tiller_delay))

        leaf = leaf_as_dict(stem_prop=stem_prop, leaf_prop=leaf_prop, rank=rank, wl=wl)
        vid_leaf = g.add_child(vid_stem, edge_type="+", **leaf)  
        tt_leaf = add_development(g=g, vid=vid_leaf, tt=tt_leaf, dtt=dtt_leaf, rate=plastochron)
        add_leaf_senescence(g=g, vid_leaf=vid_leaf, leaf_lifespan=leaf_lifespan, end_juv=end_juv)


    # Tillers

    for i in range(nb_tillers):
        # add a tiller
        vid, time = tiller_points.pop(0)

        new_tillers = add_tiller(g=g, vid=vid, start_time=time, phyllochron=phyllochron, plastochron=plastochron, tiller_delay=tiller_delay, 
                                stem_duration=stem_duration, leaf_duration=leaf_duration, leaf_lifespan=leaf_lifespan, end_juv=end_juv, 
                                reduction_factor=reduction_factor, 
                                height=height, leaf_area=leaf_area, wl=wl, diam_base=diam_base, diam_top=diam_top,
                                insertion_angle=insertion_angle, scurv=scurv, curvature=curvature,
                                klig=klig, swmax=swmax, f1=f1, f2=f2,
                                stem_q=stem_q, rmax=rmax, skew=skew, 
                                phyllotactic_angle=phyllotactic_angle, phyllotactic_deviation=phyllotactic_deviation,
                                tiller_angle=tiller_angle, gravitropism_coefficient=gravitropism_coefficient,
                                plant_orientation=plant_orientation,
                                spiral=True,
                                nb_short_phy=nb_short_phy,
                                short_phy_height=short_phy_height) 

        tiller_points.extend(new_tillers)
        # Here we consider that the list is sorted by the time
        tiller_points.sort(key=lambda x: x[1])
        tiller_points = tiller_points[:nb_tillers-i]



    g = fat_mtg(g)

    # Compute geometry
    g = mtg_interpreter(g, classic=classic)

    return g  # noqa: RET504


def majors_axes_regression(x, y):
    """Performs a major axis regression
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
    y = (a * (-b * xo + a * yo) - b * c) / (a**2 + b**2)
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
    return list(zip(*list(map(tuple, (leaf_x, leaf_y, leaf_z, leaf_r)))))


def as_leaf(polyline):
    """Compute leaf x,y,s,r tuple, length, radius max and azimuth from a
    polyline (reverse of 'as_polyline')

    Args:
        polyline: a [(x,y,z,r), ...] list of tuple sampling leaf midrib

    Returns:
        (x, y, s, r), length, radius_max , azimuth, origin of the leaf
    """

    x, y, z, r = [np.array(x).astype(float) for x in list(zip(*polyline))]
    xo, yo, zo = x[0], y[0], z[0]
    sx = curvilinear_abscisse(x, y, z)
    a, b, c = majors_axes_regression(x, y)
    
    length = sx.max()
    radius_max = r.max()
    s = sx / length
    r /= radius_max
    origin = (xo, yo, zo)
    y_leaf = z - zo
    xp, yp = list(zip(*[line_projection(a, b, c, x[0], x[1]) for x in zip(x, y)]))
    x_leaf = curvilinear_abscisse(xp, yp)
    sxp = curvilinear_abscisse(x_leaf, y_leaf)
    azimuth = np.degrees(np.arctan2(yp[-1] - yp[0], xp[-1] - xp[0]))
    return (
        (x_leaf / sxp.max(), y_leaf / sxp.max(), s, r),
        length,
        radius_max,
        azimuth,
        origin,
    )


def as_json(plant):
    """convert a plant dimension + shape table in a cereals json input
    (reverse of 'as_plant')"""
    internode = plant.L_internode.values  # noqa: PD011
    diameter = plant.W_internode.values  # noqa: PD011
    stem = [0, *internode.cumsum().tolist()]
    stem_diameter = [diameter[0]] + diameter.tolist()  # noqa: RUF005
    polylines = [
        as_polyline(leaf, length, width, (0, 0, h), azim)
        for leaf, length, width, h, azim in zip(
            plant.leaf_shape,
            plant.L_blade,
            plant.W_blade,
            plant.h_ins,
            plant.leaf_azimuth,
        )
    ]

    return {
        "leaf_polylines": polylines,
        "leaf_order": plant.leaf_rank.values.tolist(),  # noqa: PD011
        "stem": [(0, 0, z, r) for z, r in zip(stem, stem_diameter)],
    }


def as_plant(json):
    """restore plant dimension + leaves representation of a plant encoded in
    json (reverse of as_json)"""

    ranks = json["leaf_order"]
    leaves, l_leaf, w_leaf, azimuth, origin = list(
        zip(*list(map(as_leaf, json["leaf_polylines"])))
    )
    # use rank as index
    df = pd.DataFrame(  # noqa: PD901
        {
            "rank": ranks,
            "l_leaf": l_leaf,
            "w_leaf": w_leaf,
            "azimuth": azimuth,
            "hins": [ori[2] for ori in origin],
        }
    ).set_index("rank")
    df = df.sort_index()  # noqa: PD901
    leaves = {rank: leaves[i] for i, rank in enumerate(ranks)}
    x_stem, y_stem, z_stem, r_stem = list(zip(*json["stem"]))
    df["diam"] = interp1d(z_stem, r_stem, bounds_error=False, fill_value=r_stem[-1])(
        df.hins
    )
    df["ff"] = [get_form_factor(leaves[r]) for r in df.index]
    df["area"] = df.l_leaf * df.w_leaf * df.ff 
    stem = [0, *df.hins.tolist()]
    df["internode"] = np.diff(stem)
    df["ntop"] = df.index.max() - df.index + 1
    # re-index leaves with ntop
    leaves = {df.ntop[rank]: leaves[rank] for rank in df.index}
    blades = pd.DataFrame(
        {
            "L_blade": df.l_leaf,
            "S_blade": df.area,
            "W_blade": df.w_leaf,
            "ntop": df.ntop,
            "plant": 1,
            "leaf_azimuth": df.azimuth,
            "form_factor": df.ff,
        }
    )

    stem = pd.DataFrame(
        {
            "L_internode": df.internode,
            "L_sheath": 0,
            "W_internode": df.diam,
            "W_sheath": df.diam,
            "h_ins": df.hins,
            "ntop": df.ntop,
            "plant": 1,
        }
    )
    return blades, stem, leaves


"""
def build_shoot(
    nb_phy,
    height,
    leaf_area,
    wl,
    diam_base,
    diam_top,
    insertion_angle,
    scurv,
    curvature,
    # alpha,
    klig, swmax, f1, f2,
    stem_q,
    rmax,
    skew,
    phyllotactic_angle,
    phyllotactic_deviation,
    plant_orientation=45,
    spiral=True,
    nb_short_phy=4
):
    '''create a shoot

    Args:
        stem_radius: (float) the stem radius
        insertion_heights: list of each leaf insertion height
        leaf_lengths: list of each leaf length (blade length)
        leaf_areas: list of each blade area
        collar_visible: list of each collar height or True if the collar is visible and False if it is not
        leaf_shapes: list of each leaf shape, if it is not known write None
        leaf_azimuths: list of each leaf azimuth, if it is not known write None

    !!! Not used anymore

    '''

    ranks = range(1, nb_phy + 1)
    ntop = max(ranks) - np.array(ranks) + 1

    ## Stem
    # internode lengths
    #nb_young_phy = int(
    #    round((nb_phy - 1.95) / 1.84 / 1.3)
    #)  # Lejeune and Bernier formula + col

    # nb_young_phy = int(
    #     round((nb_phy - 1.95) / 1.84 / 1.3)
    # )
    # nb_short_phy = 4
    short_phy_height = 2

    pseudostem_height = nb_short_phy * short_phy_height

    pseudostem = np.array([short_phy_height*i for i in range(1, nb_short_phy+1)])
    stem = np.array(geometric_dist(height, nb_phy-nb_short_phy, q=stem_q, u0=pseudostem_height))
    insertion_heights = np.concatenate((pseudostem, stem), axis=0)
    # stem = np.array([pseudostem_height+i*(height - pseudostem_height)/(nb_phy - nb_young_phy) for i in range(nb_young_phy+1,nb_phy+1)])
    # insertion_heights = np.array(geometric_dist(height, nb_phy, q=stem_q)) #, u0=young_phy_height))
    # insertion_heights = np.array(collar_heights_kaitaniemi(height, nb_phy))

    # stem diameters
    # stem_diameters = [diam_base] * nb_young_phy + np.linspace(
    #     diam_base, diam_top, nb_phy - nb_young_phy
    # ).tolist()
    stem_diameters = np.linspace(diam_base, diam_top, nb_phy).tolist()

    stem = stem_dimension(
        h_ins=insertion_heights, d_internode=stem_diameters, ntop=ntop
    )

    ## Leaves

    # leaf length repartition along axis
    # leaf_areas_stem = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))
    # leaf_areas_pseudostem = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))
    # leaf_areas = np.concatenate((leaf_areas_pseudostem, leaf_areas_stem), axis=0)
    leaf_areas = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))


    # leaf shapes
    a_leaf = parametric_leaf(
        nb_segment=10,
        insertion_angle=insertion_angle,
        scurv=scurv,
        curvature=curvature,
        klig=klig, swmax=swmax, f1=f1, f2=f2
    )

    leaf_shapes = [a_leaf] * nb_phy 

    blades = blade_dimension(area=leaf_areas, ntop=ntop, wl=wl)

    # tck = shape_to_surface(a_leaf, wl)

    # ff = [get_form_factor(leaf) for leaf in leaf_shapes]

    # phyllotaxy
    leaf_azimuths = leaf_azimuth(
        size=nb_phy,
        phyllotactic_angle=phyllotactic_angle,
        phyllotactic_deviation=phyllotactic_deviation,
        plant_orientation=plant_orientation,
        spiral=spiral,
    )
    leaf_azimuths[1:] = np.diff(leaf_azimuths)

    ## df
    axis = blades + stem
    axis["leaf_azimuth"] = leaf_azimuths
    axis["leaf_rank"] = ranks
    axis["leaf_shape"] = [leaf_shapes[n - 1] for n in ranks]
    # df["wl"] = [wl for _ in df.leaf_rank]
    # df["leaf_tck"] = [(tck) for _ in df.leaf_rank]
    return axis, cereals_generator(plant=axis)


def shoot_at_stage(shoot, stage):
    df = shoot.loc[shoot["leaf_rank"] <= stage, :]  # noqa: PD901
    return df, cereals_generator(plant=df)
"""