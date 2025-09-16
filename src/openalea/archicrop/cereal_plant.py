"""Generate a geometric-based MTG representation of a cereal plant"""

from __future__ import annotations

import numpy as np

from openalea.mtg import MTG, fat_mtg
from openalea.mtg.turtle import TurtleFrame

from .cereal_axis import bell_shaped_dist, geometric_dist
from .cereal_leaf import parametric_cereal_leaf, shape_to_surface
from .development import add_development, compute_potential_growth_rate
from .geometry import compute_organ_geometry
from .internode import stem_as_dict, stem_dimension
from .leaf import leaf_as_dict, leaf_dimension
from .plant import leaf_azimuth
from .turtle import CerealsTurtle


def cereal_stem_properties(nb_phy, nb_short_phy, short_phy_height, height, stem_q, diam_top, diam_base, ntop):
    '''
    Return dict of stem element properties
    '''
    # short_phy_height = 2
    pseudostem_height = nb_short_phy * short_phy_height

    pseudostem = np.array([short_phy_height*i for i in range(1, nb_short_phy+1)]) if nb_short_phy > 0 else np.array([])
    stem = np.array(geometric_dist(height-pseudostem_height, nb_phy-nb_short_phy, q=stem_q, u0=short_phy_height))
    stem = np.array([h+pseudostem_height for h in stem])
    insertion_heights = np.concatenate((pseudostem, stem), axis=0) if nb_short_phy > 0 else stem

    stem_diameters = np.linspace(diam_base, diam_top, nb_phy).tolist()

    return stem_dimension(
        h_ins=insertion_heights, d_internode=stem_diameters, ntop=ntop
    )


def cereal_leaf_properties(nb_phy, leaf_area, rmax, skew, insertion_angle, scurv, curvature, 
                    klig, swmax, f1, f2, ntop, wl, phyllotactic_angle, phyllotactic_deviation, plant_orientation, spiral):
    '''
    Return dict of blade properties
    '''
    # leaf areas
    leaf_areas = np.array(bell_shaped_dist(leaf_area, nb_phy, rmax, skew))

    # leaf shapes
    a_leaf = parametric_cereal_leaf(
        nb_segment=10,
        insertion_angle=insertion_angle,
        scurv=scurv,
        curvature=curvature,
        klig=klig, swmax=swmax, f1=f1, f2=f2
    )

    blades = leaf_dimension(area=leaf_areas, ntop=ntop, wl=wl)

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


def add_leaf_senescence(g, vid_leaf, leaf_lifespan, end_juv):
    """Add leaf senescence time to the leaf vertex"""
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
    nb_short_phy = max(nb_short_phy - rank, 0)

    ranks = range(1, nb_phy + 1)
    ntop = max(ranks) - np.array(ranks) + 1

    # see how properties change with reduction factor and order p.r^o
    stem_prop = cereal_stem_properties(nb_phy, nb_short_phy, short_phy_height, height, stem_q, diam_top, diam_base, ntop)
    leaf_prop = cereal_leaf_properties(nb_phy, leaf_area, rmax, skew, insertion_angle, scurv, curvature, 
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
        compute_potential_growth_rate(g=g, vid=vid_stem)

        leaf = leaf_as_dict(stem_prop=stem_prop, leaf_prop=leaf_prop, rank=rank, wl=wl)
        vid_leaf = g.add_child(vid_stem, edge_type="+", **leaf)
        tt_leaf = add_development(g=g, vid=vid_leaf, tt=tt_leaf, dtt=dtt_leaf, rate=plastochron)
        compute_potential_growth_rate(g=g, vid=vid_leaf)
        add_leaf_senescence(g=g, vid_leaf=vid_leaf, leaf_lifespan=leaf_lifespan, end_juv=end_juv)
        
        tillers.append((vid_stem, tt_stem + tiller_delay))

    return tillers


def mtg_interpreter(g, classic=False):
    """
    Compute/update the geometry on each node of the MTG using Turtle geometry.
    """
    # BUG : sub_mtg mange le vertex plant => on perd la plante !
    # plants = g.component_roots_at_scale(g.root, scale=1)
    # nplants = g.nb_vertices(scale=1)
    # gt = MTG()

    # for plant in plants:
    #   gplant = g.sub_mtg(plant)
    turtle = CerealsTurtle()
    visitor = CerealsVisitor(classic)

    scene = TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False, all_roots=True)  # noqa: F841

    return g



class CerealsVisitor:
    """Performs geometric interpretation of mtg nodes"""

    def __init__(self, classic):
        self.classic = classic

    def __call__(self, g, v, turtle):
        # 1. retrieve the node
        n = g.node(v)

        # Go to plant position if first plant element
        if n.parent() is None:
            turtle.move(0, 0, 0)
            # initial position to be compatible with canMTG positioning
            turtle.setHead(0, 0, 1, -1, 0, 0)

        # Manage inclination of tiller
        if g.edge_type(v) == "+" and not n.label.startswith("Leaf"):
            # axis_id = g.complex(vid); g.property('insertion_angle')
            # TODO : vary as function of age and species(e.g. rice)
            # print(n.label, n.visible_length, n.tiller_angle)
            angle = 2*n.tiller_angle if g.order(v) == 1 else n.tiller_angle 
            turtle.down(angle)
            turtle.elasticity = n.gravitropism_coefficient 
            turtle.tropism = (0, 0, 1)

        # incline turtle at the base of stems,
        if n.label.startswith("Stem"):
            azim = float(n.azimuth) if n.azimuth else 0.0
            if azim:
                # print 'node', n._vid, 'azim ', azim
                turtle.rollL(azim)

        if n.label.startswith("Leaf") or n.label.startswith("Stem"):  # noqa: SIM102
            # update geometry of elements
            if n.length > 0:
                mesh = compute_organ_geometry(n, self.classic)
                if mesh:  # To DO : reset to None if calculated so ?
                    n.geometry = turtle.transform(mesh)
                    n.anchor_point = turtle.getPosition()
                    n.heading = turtle.getHeading()

        # 3. Update the turtle and context
        turtle.setId(v)
        if n.label.startswith("Stem"):
            if n.visible_length > 0:  
                turtle.f(n.visible_length)
            turtle.context.update({"top": turtle.getFrame()})
        if n.label.startswith("Leaf"):  # noqa: SIM102
            if n.lrolled > 0:
                turtle.f(n.lrolled)
                turtle.context.update({"top": turtle.getFrame()})


def cereal(nb_phy, phyllochron, plastochron, stem_duration, leaf_duration, 
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
    stem_prop = cereal_stem_properties(nb_phy=nb_phy, nb_short_phy=nb_short_phy, short_phy_height=0.01, height=1, 
                                        stem_q=stem_q, diam_top=diam_top, diam_base=diam_base, ntop=ntop)
    leaf_prop = cereal_leaf_properties(nb_phy=nb_phy, leaf_area=1, rmax=rmax, skew=skew, insertion_angle=insertion_angle, 
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
        # compute_potential_growth_rate(g=g, vid=vid_stem)
        tiller_points.append((vid_stem, tt_stem + tiller_delay))

        leaf = leaf_as_dict(stem_prop=stem_prop, leaf_prop=leaf_prop, rank=rank, wl=wl)
        vid_leaf = g.add_child(vid_stem, edge_type="+", **leaf)  
        tt_leaf = add_development(g=g, vid=vid_leaf, tt=tt_leaf, dtt=dtt_leaf, rate=plastochron)
        # compute_potential_growth_rate(g=g, vid=vid_leaf)
        add_leaf_senescence(g=g, vid_leaf=vid_leaf, leaf_lifespan=leaf_lifespan, end_juv=end_juv)


    # Tillers

    for i in range(nb_tillers):
        # add a tiller
        vid, time = tiller_points.pop(0)

        new_tillers = add_tiller(g=g, vid=vid, start_time=time, phyllochron=phyllochron, plastochron=plastochron, tiller_delay=tiller_delay, 
                                stem_duration=stem_duration, leaf_duration=leaf_duration, leaf_lifespan=leaf_lifespan, end_juv=end_juv, 
                                reduction_factor=reduction_factor, 
                                height=1, leaf_area=1, wl=wl, diam_base=diam_base, diam_top=diam_top,
                                insertion_angle=insertion_angle, scurv=scurv, curvature=curvature,
                                klig=klig, swmax=swmax, f1=f1, f2=f2,
                                stem_q=stem_q, rmax=rmax, skew=skew, 
                                phyllotactic_angle=phyllotactic_angle, phyllotactic_deviation=phyllotactic_deviation,
                                tiller_angle=tiller_angle, gravitropism_coefficient=gravitropism_coefficient,
                                plant_orientation=plant_orientation,
                                spiral=True,
                                nb_short_phy=nb_short_phy,
                                short_phy_height=0.01) 

        tiller_points.extend(new_tillers)
        # Here we consider that the list is sorted by the time
        tiller_points.sort(key=lambda x: x[1])
        tiller_points = tiller_points[:nb_tillers-i]


    # Find number of phytomers and order of each axis 
    # Determine a ratio of total leaf area to allocate to each axis (with or without reduction factor)
    # Browse all leaves of all axes and multiply them by total_leaf_area * ratio_axis

    # Find number of phytomers and order of each axis 
    axes = g.vertices(scale=2)  # scale=2: axes
    main_axis = axes[0]
    axis_orders = {axis: g.order(axis) for axis in axes}

    # Compute the number of leaves per axis
    leaves_per_axis = {
        axis: 
        [
            vid for vid in g.components(axis)
            if g.node(vid).label.startswith("Leaf")
        ]
        for axis in axes
    }
    nb_leaves_per_axis = {axis: len(leaves) for axis, leaves in leaves_per_axis.items()}

    stem_elements_per_axis = {
        axis: 
        [
            vid for vid in g.components(axis)
            if g.node(vid).label.startswith("Stem")
        ]
        for axis in axes
    }
    # total_nb_leaves = sum(nb_leaves_per_axis.values())

    # Compute reduction factor for each axis (main axis: 1, tillers: reduction_factor^order)
    axis_factors = {
        axis: (1.0 * nb_leaves_per_axis[axis]
               if axis == main_axis 
               else (reduction_factor ** axis_orders[axis]) * nb_leaves_per_axis[axis])  
        for axis in axes
    }
    total_factor = sum([axis_factors[axis] for axis in axes])

    # Allocate total leaf area to each axis proportionally
    for axis in axes:
        leaves = leaves_per_axis[axis]
        stem = stem_elements_per_axis[axis]
        factor = axis_factors[axis]
        # Area per leaf for this axis
        area_per_leaf = leaf_area * factor / total_factor # if nb_leaves_per_axis[axis] > 0 else 0
        height_per_axis = height * factor / total_factor # if nb_leaves_per_axis[axis] > 0 else 0
        for vid in leaves:
            scaled_leaf_area = g.node(vid).leaf_area * area_per_leaf
            g.node(vid).leaf_area = scaled_leaf_area  
            # g.node(vid).visible_leaf_area *= area_per_leaf
            g.node(vid).mature_length = np.sqrt(scaled_leaf_area / 0.75 / wl)
            g.node(vid).length = np.sqrt(scaled_leaf_area / 0.75 / wl) 
            g.node(vid).visible_length = np.sqrt(scaled_leaf_area / 0.75 / wl) 
            g.node(vid).shape_max_width = g.node(vid).length * wl
            compute_potential_growth_rate(g=g, vid=vid)
        for vid in stem:
            scaled_stem_length = g.node(vid).mature_length * height_per_axis
            g.node(vid).mature_length = scaled_stem_length
            g.node(vid).length = scaled_stem_length
            g.node(vid).visible_length = scaled_stem_length
            compute_potential_growth_rate(g=g, vid=vid)


    g = fat_mtg(g)

    # Compute geometry
    g = mtg_interpreter(g, classic=classic)

    return g  # noqa: RET504

