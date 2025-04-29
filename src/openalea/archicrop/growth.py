from __future__ import annotations

from scipy.interpolate import splev

from openalea.mtg.traversal import pre_order2_with_filter

from .geometry import (
    CerealsTurtle,
    CerealsVisitor,
    TurtleFrame,
    addSets,
    leaf_mesh,
    stem_mesh,
)


def init_visible_variables(g):
    """Initialize leaf and stem visible lengths"""

    for k in g.properties()["visible_length"]:
        g.properties()["visible_length"][k] = 0.0

    return g


def equal_dist(increment, growing_organs):
    '''return initial distribution per organ (H: same amount for all organs)'''
    return {vid: increment/len(growing_organs) for vid in growing_organs}


def demand_dist(increment, growing_organs):
    '''return distribution of increment per organ proportionnal to potential length or area'''
    sum_growing_organs = sum([value["potential"] for value in growing_organs.values()])
    return {vid: increment*values["potential"]/sum_growing_organs 
            for vid, values in growing_organs.items()}


def get_growing_and_senescing_organs_potential_visible(g, time, prev_time):
    """Identify growing organs and their potential at a given time"""

    growing_internodes = {}
    growing_leaves = {}
    senescing_leaves = {}

    for vid, la in g.properties()["leaf_area"].items(): 
        n = g.node(vid)
        if n.start_tt <= time <= n.end_tt or prev_time < n.end_tt < time:
            growing_leaves[vid] = {"potential": la, "visible": n.visible_leaf_area}
        elif n.senescence <= time and not n.dead: # and n.srt > 0: 
            senescing_leaves[vid] = {"potential": n.visible_leaf_area, "visible": n.senescent_area}

    for vid, ml in g.properties()["mature_length"].items(): 
        n = g.node(vid)
        if n.label.startswith("Stem") and (n.start_tt <= time < n.end_tt or prev_time < n.end_tt < time):
            growing_internodes[vid] = {"potential": ml, "visible": n.visible_length}

    return growing_internodes, growing_leaves, senescing_leaves


def get_growing_organs(g, time, prev_time):
    """Identify growing organs and their potential at a given time"""

    growing_internodes = []
    growing_leaves = []

    for vid in g.properties()["leaf_area"]: 
        n = g.node(vid)
        if n.start_tt <= time <= n.end_tt or prev_time < n.end_tt < time:
            growing_leaves.append(vid)

    for vid in g.properties()["mature_length"]: 
        n = g.node(vid)
        if n.label.startswith("Stem") and (n.start_tt <= time < n.end_tt or prev_time < n.end_tt < time):
            growing_internodes.append(vid)

    return growing_internodes, growing_leaves


def get_senescing_organs(g, time):
    """Identify senescing blades at a given time"""

    senescing_leaves = []
    for vid in g.properties()["leaf_area"]: 
        n = g.node(vid)
        if n.senescence <= time and n.srt > 0: 
            senescing_leaves.append(vid)
    return senescing_leaves


def distribute_to_potential(growing_organs, increment_to_distribute, distribution_function):
    "Distribute increment among growing organs up to potential of each organ"

    increment_for_each_organ = dict.fromkeys(growing_organs, 0.0)

    while len(growing_organs) > 0 and increment_to_distribute > 1e-4: 
        incr_temp = distribution_function(increment_to_distribute, growing_organs)
        increment_to_distribute = 0.0

        for vid, val in incr_temp.items():
        #     increment_for_each_organ[vid] += val # added here
        # for vid in increment_for_each_organ.keys():
            if vid in growing_organs:
                potential_increment = min(val, growing_organs[vid]["potential"] - growing_organs[vid]["visible"])
                increment_for_each_organ[vid] += potential_increment
                growing_organs[vid]["visible"] += potential_increment
                increment_to_distribute += val - potential_increment
                # if incr_temp[vid] + growing_organs[vid]["visible"] > growing_organs[vid]["potential"]: # not good
                    # remaining = increment_for_each_organ[vid] + growing_organs[vid]["visible"] - growing_organs[vid]["potential"]
                    # increment_to_distribute += remaining
                    # increment_for_each_organ[vid] -= remaining
                # growing_organs[vid]["visible"] += increment_for_each_organ[vid] # and added here

        # filter growing_organs
        growing_organs = {vid: {"potential": values["potential"], "visible": values["visible"]} 
                          for vid, values in growing_organs.items() 
                          if values["visible"] < values["potential"]}

    return increment_for_each_organ



def distribute_among_organs(g, time, prev_time, daily_dynamics):
    """Fill the dictionnary of variables recording growth process at a given time"""

    # compute number of growing organs
    growing_internodes, growing_leaves, senescing_leaves = get_growing_and_senescing_organs_potential_visible(g, time, prev_time)

    # compute number of senscent organs
    # senescing_leaves = get_senescing_organs(g, time)

    # set amount of growth to distribute among organs
    height_to_distribute = daily_dynamics["Height increment"] 
    LA_to_distribute = daily_dynamics["Leaf area increment"] 
    sen_LA_to_distribute = daily_dynamics["Senescent leaf area increment"]
    
    height_for_each_internode = dict.fromkeys(growing_internodes, 0.0)
    LA_for_each_leaf = dict.fromkeys(growing_leaves, 0.0)
    sen_LA_for_each_leaf = dict.fromkeys(senescing_leaves, 0.0)


    if LA_to_distribute > 0.0:
        LA_for_each_leaf = distribute_to_potential(growing_leaves, LA_to_distribute, demand_dist)
    if height_to_distribute > 0.0:
        height_for_each_internode = distribute_to_potential(growing_internodes, height_to_distribute, demand_dist)
    if sen_LA_to_distribute > 0.0:
        sen_LA_for_each_leaf = distribute_to_potential(senescing_leaves, sen_LA_to_distribute, equal_dist)

    return {"height_to_distribute": height_to_distribute,
              "height_for_each_internode": height_for_each_internode,
              "LA_to_distribute": LA_to_distribute,
              "LA_for_each_leaf": LA_for_each_leaf,
              "sen_LA_to_distribute": sen_LA_to_distribute,
              "sen_LA_for_each_leaf": sen_LA_for_each_leaf,
              "growing_internodes": growing_internodes,
              "growing_leaves": growing_leaves,
              "senescing_leaves": senescing_leaves}


def compute_organ(vid, element_node, time, growth, classic=False):
    """Compute geometry of Adel base elements (LeafElement and StemElement).
    element_node should be a mtg node proxy."""

    n = element_node
    geom = None

    if n.label.startswith("Leaf") or n.label.startswith("Stem"):
        if time < n.start_tt:
            n.visible_length = 0.0
            n.grow = False
            if n.label.startswith("Leaf"):
                n.leaf_lengths.append(0.0)
                n.senescent_lengths.append(0.0)
            elif n.label.startswith("Stem"):
                n.stem_lengths.append(0.0)
        elif n.start_tt <= time:
            n.age = time - n.start_tt

            if time <= n.end_tt:
                n.grow = True
                if n.visible_length < n.mature_length:
                    if n.label.startswith("Leaf") and growth["LA_to_distribute"] > 0.0:
                        LA_for_this_leaf = growth["LA_for_each_leaf"][vid]
                        update_leaf_growth(n, LA_for_this_leaf)
                    elif n.label.startswith("Stem") and growth["height_to_distribute"] > 0.0:
                        height_for_this_internode = growth["height_for_each_internode"][vid]
                        update_stem_growth(n, height_for_this_internode)
            elif time > n.end_tt:
                if n.label.startswith("Leaf") and growth["sen_LA_to_distribute"] > 0.0:
                    if n.senescence <= time and not n.dead and vid in growth["sen_LA_for_each_leaf"]:
                        sen_LA_for_this_leaf = growth["sen_LA_for_each_leaf"][vid]
                        update_leaf_senescence(n, sen_LA_for_this_leaf)
                    else:
                        n.senescent_lengths.append(n.visible_length)
                else:
                    finalize_organ_growth(n)

        # Evolution of organ geometry with age
        stem_diameter = min(n.stem_diameter/2 * (1+0.5*(time - n.start_tt) / (n.end_tt - n.start_tt)), n.stem_diameter)
        leaf_inclination = min(1.5*(time - n.start_tt) / (n.end_tt - n.start_tt), 2)

    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
            if n.shape is not None and n.srb is not None:
                geom = leaf_mesh(
                    leaf=n.shape,
                    L_shape=n.mature_length,
                    Lw_shape=n.shape_max_width,
                    length=n.visible_length,
                    s_base=n.srb,
                    s_top=n.srt,
                    # flipx allows x-> -x to place the shape along
                    # with the tiller positioned with
                    # turtle.down()
                    flipx=True,
                    inclination=leaf_inclination, 
                    stem_diameter=stem_diameter,
                )

            ### Add a mesh that is senescent
            if n.senescent_length > 0.0001: # filter less than 0.001 mm leaves  # noqa: SIM102
                if n.srt is not None: # and n.shape is not None
                    geom_senescent = leaf_mesh(
                        leaf=n.shape,
                        L_shape=n.mature_length,
                        Lw_shape=n.shape_max_width,
                        length=n.visible_length,
                        s_base=n.srt,
                        s_top=1.,
                        # flipx allows x-> -x to place the shape along
                        # with the tiller positioned with
                        # turtle.down()
                        flipx=True,
                        inclination=leaf_inclination, 
                        stem_diameter=stem_diameter,
                    )
                    n.geometry_senescent = geom_senescent


            if n.lrolled > 0:
                rolled = stem_mesh(
                    n.lrolled, n.lrolled, n.d_rolled, classic
                )
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        geom = stem_mesh(n.length, n.visible_length, stem_diameter, stem_diameter, classic)


    return geom, growth

def update_leaf_growth(n, LA_for_this_leaf):
    if n.visible_leaf_area + LA_for_this_leaf < n.leaf_area:
        n.visible_leaf_area += LA_for_this_leaf
        relative_visible_area = n.visible_leaf_area / n.leaf_area
        # tck = shape_to_surface(n.shape, n.wl)
        n.visible_length = float(splev(x=relative_visible_area, tck=n.tck)) * n.mature_length
    else:
        n.visible_leaf_area = n.leaf_area
        n.visible_length = n.mature_length
    n.leaf_lengths.append(n.visible_length)
    n.senescent_lengths.append(0.0)

def update_stem_growth(n, height_for_this_internode):
    if n.visible_length + height_for_this_internode <= n.mature_length:
        n.visible_length += height_for_this_internode
    else:
        n.visible_length = n.mature_length
    n.stem_lengths.append(n.visible_length)

def update_leaf_senescence(n, sen_LA_for_this_leaf):
    n.grow = False
    if n.senescent_area + sen_LA_for_this_leaf < n.visible_leaf_area:
        n.senescent_area += sen_LA_for_this_leaf
        relative_senescent_area = n.senescent_area / n.visible_leaf_area
        # tck = shape_to_surface(n.shape, n.wl)
        n.srt = 1 - float(splev(x=relative_senescent_area, tck=n.tck))
        n.senescent_length = (1 - n.srt) * n.visible_length
    else:
        n.senescent_area = n.visible_leaf_area
        n.senescent_length = n.visible_length
        n.dead = True
    n.senescent_lengths.append(n.senescent_length)

def finalize_organ_growth(n):
    # if n.visible_length < n.mature_length:
    #     n.visible_length = n.visible_length
    # else:
    #     n.visible_length = n.mature_length
    if n.label.startswith("Leaf"):
        n.leaf_lengths.append(n.visible_length)
        n.senescent_lengths.append(n.senescent_length)
    elif n.label.startswith("Stem"):
        n.stem_lengths.append(n.visible_length)

        

class CerealsVisitorConstrained(CerealsVisitor):
    def __init__(self, classic):
        super().__init__(classic)

    def __call__(self, g, v, turtle, time, growth):
        # 1. retrieve the node
        geoms_senesc = g.property("geometry_senescent")
        n = g.node(v)

        # Go to plant position if first plant element
        if n.parent() is None:
            turtle.move(0, 0, 0)
            # initial position to be compatible with canMTG positioning
            turtle.setHead(0, 0, 1, -1, 0, 0)

        # incline turtle at the base of stems,
        if n.label.startswith("Stem"):
            azim = float(n.azimuth) if n.azimuth else 0.0
            if azim:
                # print 'node', n._vid, 'azim ', azim
                turtle.rollL(azim)

        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            # update geometry of elements
            # if n.length > 0:
            # print(v)
            # nb_of_growing_internodes, nb_of_growing_leaves = growing_organs(g, time)
            mesh, growth = compute_organ(v, n, time, growth, self.classic)
            if mesh:  # To DO : reset to None if calculated so ?
                n.geometry = turtle.transform(mesh)
                n.anchor_point = turtle.getPosition()
            if v in geoms_senesc and geoms_senesc[v] is not None:
                geoms_senesc[v] = turtle.transform(geoms_senesc[v])

        # 3. Update the turtle and context
        turtle.setId(v)
        if n.label.startswith("Stem"):
            if n.length > 0:
                turtle.f(n.visible_length)
            turtle.context.update({"top": turtle.getFrame()})
        if n.label.startswith("Leaf"):  # noqa: SIM102
            if n.lrolled > 0:
                turtle.f(n.lrolled)
                turtle.context.update({"top": turtle.getFrame()})

        return growth
    

        


def mtg_turtle_time_with_constraint(g, time, prev_time, daily_dynamics, update_visitor=None):  # noqa: ARG001
    """Compute the geometry on each node of the MTG using Turtle geometry.

    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.

    """

    g.properties()["geometry"] = {}
    g.properties()["_plant_translation"] = {}

    max_scale = g.max_scale()

    growth = distribute_among_organs(g, time, prev_time, daily_dynamics)

    cereal_visitor = CerealsVisitorConstrained(False)

    # for id in g.vertices():
    #     print(g[id])
    
    
    def traverse_with_turtle_time(g, vid, time, growth, visitor=cereal_visitor):
        turtle = CerealsTurtle()

        def push_turtle(v):
            n = g.node(v)
            # if 'Leaf' in n.label:
            #    return False
            try:
                start_tt = n.start_tt
                if start_tt > time:
                    return False
            except:  # noqa: E722
                pass
            if g.edge_type(v) == "+":
                turtle.push()
            return True

        def pop_turtle(v):
            n = g.node(v)
            try:
                start_tt = n.start_tt
                if start_tt > time:
                    return False
            except:  # noqa: E722
                pass
            if g.edge_type(v) == "+":  # noqa: RET503
                turtle.pop()  # noqa: RET503

        if g.node(vid).label.startswith("Leaf") or g.node(vid).label.startswith("Stem"):  # noqa: SIM102
            if g.node(vid).start_tt <= time:
                growth = visitor(g, vid, turtle, time, growth)
                # turtle.push()
        # plant_id = g.complex_at_scale(vid, scale=1)

        for v in pre_order2_with_filter(g, vid, None, push_turtle, pop_turtle):
            if v == vid:
                continue
            # Done for the leaves
            if g.node(v).start_tt > time:
                continue
            growth = visitor(g, v, turtle, time, growth)

        # scene = turtle.getScene()
        return g


    for plant_id in g.component_roots_at_scale_iter(g.root, scale=max_scale):
        g = traverse_with_turtle_time(g, plant_id, time, growth)
    return g


def mtg_interpreter_with_constraint(g, classic=False):
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
    visitor = CerealsVisitorConstrained(classic)

    scene = TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False, all_roots=True)  # noqa: F841

    return g



