from __future__ import annotations

import math
import numpy as np
import sympy as sp
from scipy.interpolate import splev
from openalea.mtg.traversal import pre_order2, pre_order2_with_filter

from .geometry import CerealsTurtle, CerealsVisitor, addSets, TurtleFrame, leaf_mesh_for_growth, stem_mesh
from .plant_shape import shape_to_surface




def thermal_time(g, phyllochron, plastochron, leaf_duration, stem_duration, leaf_lifespan, end_juv): # durations 1.6
    """
    Add dynamic properties on the mtg to simulate development

    :param g: MTG, MTG of a plant
    :param phyllochron: float, phyllochron, i.e. internode appearance rate (in °C.day/internode)
    :param plastochron: float, plastochron, i.e. leaf appearance rate (in °C.day/leaf)
    :param leaf_duration: float, phyllochronic time for a leaf to develop from tip appearance to collar appearance (/phyllochron)
    :param stem_duration: float, phyllochronic time for a stem to develop from base to top (/phyllochron)


    :return: MTG
    """

    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        tt_stem = 0
        tt_leaf = 0
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))
        # stem_ids = g.Trunk(v)
        # nb_stems = len(stem_ids)
        # nb_sectors = 1
        dtt_stem = phyllochron * stem_duration
        dtt_leaf = plastochron * leaf_duration

        for metamer in pre_order2(g, v):
            nm = g.node(metamer)

            if "Stem" in nm.label:
                nm.start_tt = tt_stem
                nm.end_tt = tt_stem + dtt_stem
                tt_stem += phyllochron
            elif "Leaf" in nm.label:
                nm.start_tt = tt_leaf
                nm.end_tt = tt_leaf + dtt_leaf
                tt_leaf += plastochron

                # check in which pheno stage is start_tt to know which value of lifespan to use
                if nm.start_tt < end_juv:
                    nm.senescence = nm.start_tt + leaf_lifespan[0] 
                else:
                    nm.senescence = nm.start_tt + leaf_lifespan[1]

    return g


def init_visible_variables(g):
    """Initialize leaf and stem visible lengths"""

    for k in g.properties()["visible_length"].keys():
        g.properties()["visible_length"][k] = 0.0

    return g


class Growth:

    def __init__(self, height_to_distribute: float, height_for_each_internode: dict[int,float], LA_to_distribute: float, LA_for_each_leaf: dict[int,float], 
                 sen_LA_to_distribute: float, sen_LA_for_each_leaf: dict[int,float]) -> None:
        self.height_to_distribute = height_to_distribute
        self.height_for_each_internode = height_for_each_internode
        self.LA_to_distribute = LA_to_distribute
        self.LA_for_each_leaf = LA_for_each_leaf
        self.sen_LA_to_distribute = sen_LA_to_distribute
        self.sen_LA_for_each_leaf = sen_LA_for_each_leaf


def equal_dist(increment, growing_organs):
    '''return initial distribution per organ (H: same amount for all organs)'''
    return {vid: increment/len(growing_organs) for vid in growing_organs.keys()}


def demand_dist(increment, growing_organs):
    '''return distribution of increment per organ proportionnal to potential length or area'''
    sum_growing_organs = sum([value["potential"] for value in growing_organs.values()])
    return {vid: increment*values["potential"]/sum_growing_organs 
            for vid, values in growing_organs.items()}


def get_growing_and_senescing_organs_potential_visible(g, time):
    """Identify growing organs and their potential at a given time"""

    growing_internodes = {}
    growing_leaves = {}
    senescing_leaves = {}

    for vid, la in g.properties()["leaf_area"].items(): 
        n = g.node(vid)
        if n.start_tt <= time < n.end_tt:
            growing_leaves[vid] = {"potential": la, "visible": n.visible_leaf_area}
        elif n.senescence <= time and not n.dead: # and n.srt > 0: 
            senescing_leaves[vid] = {"potential": n.visible_leaf_area, "visible": n.senescent_area}

    for vid, ml in g.properties()["mature_length"].items(): 
        n = g.node(vid)
        if n.label.startswith("Stem") and n.start_tt <= time < n.end_tt:
            growing_internodes[vid] = {"potential": ml, "visible": n.visible_length}

    return growing_internodes, growing_leaves, senescing_leaves


def get_growing_organs(g, time):
    """Identify growing organs and their potential at a given time"""

    growing_internodes = []
    growing_leaves = []

    for vid in g.properties()["leaf_area"].keys(): 
        n = g.node(vid)
        if n.start_tt <= time < n.end_tt:
            growing_leaves.append(vid)

    for vid in g.properties()["mature_length"].keys(): 
        n = g.node(vid)
        if n.label.startswith("Stem") and n.start_tt <= time < n.end_tt:
            growing_internodes.append(vid)

    return growing_internodes, growing_leaves


def get_senescing_organs(g, time):
    """Identify senescing blades at a given time"""

    senescing_leaves = []
    for vid in g.properties()["leaf_area"].keys(): 
        n = g.node(vid)
        if n.senescence <= time and n.srt > 0: 
            senescing_leaves.append(vid)
    return senescing_leaves


def distribute_to_potential(growing_organs, increment_to_distribute, distribution_function):
    "Distribute increment among growing organs up to potential of each organ"

    increment_for_each_organ = {vid: 0.0 for vid in growing_organs}

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



def distribute_among_organs(g, time, increments):
    """Fill the dictionnary of variables recording growth process at a given time"""

    # compute number of growing organs
    growing_internodes, growing_leaves, senescing_leaves = get_growing_and_senescing_organs_potential_visible(g, time)

    # compute number of senscent organs
    # senescing_leaves = get_senescing_organs(g, time)

    # set amount of growth to distribute among organs
    height_to_distribute = increments["Height increment"] 
    LA_to_distribute = increments["Leaf area increment"] 
    sen_LA_to_distribute = increments["Senescent leaf area increment"]
    
    height_for_each_internode = {vid: 0.0 for vid in growing_internodes}
    LA_for_each_leaf = {vid: 0.0 for vid in growing_leaves}
    sen_LA_for_each_leaf = {vid: 0.0 for vid in senescing_leaves}

    # print(time)
    # print(sen_LA_to_distribute)

    if LA_to_distribute > 0.0:
        LA_for_each_leaf = distribute_to_potential(growing_leaves, LA_to_distribute, demand_dist)
        # print(LA_for_each_leaf)
    if height_to_distribute > 0.0:
        height_for_each_internode = distribute_to_potential(growing_internodes, height_to_distribute, demand_dist)
    if sen_LA_to_distribute > 0.0:
        sen_LA_for_each_leaf = distribute_to_potential(senescing_leaves, sen_LA_to_distribute, equal_dist)
        # print(time)
    #     print(sen_LA_to_distribute == sum(list(sen_LA_for_each_leaf.values())))
    # print('')

    growth = {"height_to_distribute": height_to_distribute,
              "height_for_each_internode": height_for_each_internode,
              "LA_to_distribute": LA_to_distribute,
              "LA_for_each_leaf": LA_for_each_leaf,
              "sen_LA_to_distribute": sen_LA_to_distribute,
              "sen_LA_for_each_leaf": sen_LA_for_each_leaf,
              "growing_internodes": growing_internodes,
              "growing_leaves": growing_leaves,
              "senescing_leaves": senescing_leaves}

    return growth


'''
def compute_organ(
    vid, element_node, time, growth, classic=False
):  # see maybe with *kwarg, **kwds, etc. for time
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""
    n = element_node
    geom = None


    if n.label.startswith("Leaf") or n.label.startswith("Stem"):

        # organ not yet appeared --> never the case because only apeared organs go through this function
        if time < n.start_tt:  
            n.visible_length = 0.0
            n.grow = False
            if n.label.startswith("Leaf"):
                n.leaf_lengths.append(0.0)
                n.senescent_lengths.append(0.0)
            elif n.label.startswith("Stem"):
                n.stem_lengths.append(0.0)
    
        # organ appeared
        elif n.start_tt <= time:
            n.grow = True
            n.age = time - n.start_tt # update age of organ
        
            # starting from the lowest (oldest) growing organ
            if time <= n.end_tt: # if n.start_tt <= time < n.end_tt: # organ growing

                if n.visible_length < n.mature_length:  

                    # for growing blades
                    if n.label.startswith("Leaf") and growth["LA_to_distribute"] > 0.0:
                        LA_for_this_leaf = growth["LA_for_each_leaf"][vid]
                        
                        if n.visible_leaf_area + LA_for_this_leaf < n.leaf_area:
                            n.visible_leaf_area += LA_for_this_leaf
                            # convert constraint on area for each leaf(t) to leaf length
                            relative_visible_area = n.visible_leaf_area / n.leaf_area
                            tck = shape_to_surface(n.shape, n.wl) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            n.visible_length = float(splev(x=relative_visible_area, tck=tck)) * n.mature_length
                        else: 
                            n.visible_leaf_area = n.leaf_area
                            n.visible_length = n.mature_length

                        n.leaf_lengths.append(n.visible_length)
                        n.senescent_lengths.append(0.0)
                        
                    # for growing stem elements
                    elif n.label.startswith("Stem") and growth["height_to_distribute"] > 0.0:
                        height_for_this_internode = growth["height_for_each_internode"][vid]

                        if n.visible_length + height_for_this_internode <= n.mature_length:
                            n.visible_length += height_for_this_internode
                        else: 
                            n.visible_length = n.mature_length

                        n.stem_lengths.append(n.visible_length)

            
            # organ reaches maturity according to defined development
            elif time > n.end_tt:

                # blade senescence
                if n.label.startswith("Leaf") and growth["sen_LA_to_distribute"] > 0.0:
                    if n.senescence <= time:
                        if not n.dead:
                            n.grow = False 
                            sen_LA_for_this_leaf = growth["sen_LA_for_each_leaf"][vid]
                            
                            if n.senescent_area + sen_LA_for_this_leaf < n.leaf_area:
                                n.senescent_area += sen_LA_for_this_leaf
                                relative_senescent_area = n.senescent_area / n.leaf_area
                                tck = shape_to_surface(n.shape, n.wl) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                n.srt = 1 - float(splev(x=relative_senescent_area, tck=tck))
                                n.senescent_length = (1 - n.srt) * n.mature_length
                            else: 
                                n.senescent_area = n.leaf_area
                                n.senescent_length = n.mature_length
                                n.dead = True

                            n.senescent_lengths.append(n.senescent_length)

                        else:
                            n.senescent_lengths.append(n.mature_length)

                else:
                    if n.visible_length < n.mature_length:  # organ reaches maturity below potential 
                        n.visible_length = n.visible_length
                    else:  # organ reaches maturity at potential
                        n.visible_length = n.mature_length
                    if n.label.startswith("Leaf"):
                        n.leaf_lengths.append(n.visible_length)
                        n.senescent_lengths.append(n.senescent_length)
                    elif n.label.startswith("Stem"):
                        n.stem_lengths.append(n.visible_length)
    
    
    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
            # TODO : test if senesc or not
            # if senesc : create a mesh in green and a mesh with senesc
            if n.shape is not None and n.srb is not None:
                geom = leaf_mesh_for_growth(
                    n.shape,
                    n.mature_length,
                    n.shape_max_width,
                    n.visible_length,
                    n.srb,
                    n.srt,
                    # flipx allows x-> -x to place the shape along
                    #  with the tiller positioned with
                    # turtle.down()
                    flipx=True,
                    inclination=min(0.5 + 0.5*(time - n.start_tt) / (n.end_tt - n.start_tt),1.5), # !!!!!!!!!!!!!!!!!
                    stem_diameter=n.stem_diameter,
                )
            if n.lrolled > 0:
                rolled = stem_mesh(
                    n.lrolled, n.lrolled, n.d_rolled, classic
                )
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        geom = stem_mesh(n.length, n.visible_length, n.stem_diameter, n.stem_diameter, classic)


    return geom, growth
'''

#########################################

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
                    if n.senescence <= time and not n.dead:
                        sen_LA_for_this_leaf = growth["sen_LA_for_each_leaf"][vid]
                        update_leaf_senescence(n, sen_LA_for_this_leaf)
                    else:
                        n.senescent_lengths.append(n.visible_length)
                else:
                    finalize_organ_growth(n)

    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
            # TODO : test if senesc or not
            # if senesc : create a mesh in green and a mesh with senesc
            if n.shape is not None and n.srb is not None:
                geom = leaf_mesh_for_growth(
                    n.shape,
                    n.mature_length,
                    n.shape_max_width,
                    n.visible_length,
                    n.srb,
                    n.srt,
                    # flipx allows x-> -x to place the shape along
                    #  with the tiller positioned with
                    # turtle.down()
                    flipx=True,
                    inclination=min(0.5 + 0.5*(time - n.start_tt) / (n.end_tt - n.start_tt),1.5), # !!!!!!!!!!!!!!!!!
                    stem_diameter=n.stem_diameter,
                )

            ### Add a mesh that is senescent
            if n.senescent_length > 0.0001:
                if n.shape is not None and n.srb is not None:
                    geom_senescent = leaf_mesh_for_growth(
                        n.shape,
                        n.mature_length,
                        n.shape_max_width,
                        n.visible_length,
                        n.srt,
                        1.,
                        # flipx allows x-> -x to place the shape along
                        #  with the tiller positioned with
                        # turtle.down()
                        flipx=True,
                        inclination=min(0.5 + 0.5*(time - n.start_tt) / (n.end_tt - n.start_tt),1.5), # !!!!!!!!!!!!!!!!!
                        stem_diameter=n.stem_diameter,
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
        geom = stem_mesh(n.length, n.visible_length, n.stem_diameter, n.stem_diameter, classic)


    return geom, growth

def update_leaf_growth(n, LA_for_this_leaf):
    if n.visible_leaf_area + LA_for_this_leaf < n.leaf_area:
        n.visible_leaf_area += LA_for_this_leaf
        relative_visible_area = n.visible_leaf_area / n.leaf_area
        tck = shape_to_surface(n.shape, n.wl)
        n.visible_length = float(splev(x=relative_visible_area, tck=tck)) * n.mature_length
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
        tck = shape_to_surface(n.shape, n.wl)
        n.srt = 1 - float(splev(x=relative_senescent_area, tck=tck))
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


###########################################
        

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
            if v in geoms_senesc:
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
    

        


def mtg_turtle_time_with_constraint(g, time, increments, update_visitor=None):
    """Compute the geometry on each node of the MTG using Turtle geometry.

    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.

    """

    g.properties()["geometry"] = {}
    g.properties()["_plant_translation"] = {}

    max_scale = g.max_scale()

    growth = distribute_among_organs(g, time, increments)

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
            except:
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
            except:
                pass
            if g.edge_type(v) == "+":  # noqa: RET503
                turtle.pop()

        if g.node(vid).label.startswith("Leaf") or g.node(vid).label.startswith("Stem"):
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

    scene = TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False, all_roots=True)

    return g




# def leaf_area_function_of_length(l, L, wl, alpha=-2.3):
#     '''returns the current leaf area (i.e. twice the integral of leaf shape) given the current length s
    
#         function obtained this way : 

#         # Define the variable and parameters
#         s = sympy.symbols('l')
#         L, wl, alpha = sympy.symbols('L, wl, alpha') 

#         # Define the scaled leaf shape function with parameters
#         # alpha = -2.3
#         beta = -2 * (alpha + sqrt(-alpha))
#         gamma = 2 * sqrt(-alpha) + alpha
#         r = wl * L (alpha * s**2 + beta * s + gamma)

#         # Reflect the function to the axis x=0.5 
#         reflected_r = r.subs(s, -s+L).simplify()

#         # Find the indefinite integral (primitive)
#         primitive_f = sympy.integrate(reflected_r, s)

#         print(2*primitive_f)
    
#     '''

#     return 2*l**2*wl*math.sqrt(-alpha) + 2*alpha*l**3*wl/(3*L)


# def compute_leaf_length_increment(LA_for_each_leaf, leaf):
#     """returns the leaf length given a leaf area"""

#     L = leaf.mature_length
#     current_leaf_length = leaf.visible_length
#     wl = leaf.shape_max_width / leaf.mature_length
    
#     s = sp.symbols('s')
    
#     # Define the equation
#     equation = LA_for_each_leaf + leaf_area_function_of_length(current_leaf_length, L, wl) - leaf_area_function_of_length(s, L, wl)
#     # print ("LA_for_each_leaf", LA_for_each_leaf)

#     # Solve the equation, restricting the domain to real numbers
#     real_solutions = sp.solveset(equation, s, domain=sp.S.Reals)
#     # print("Sol", real_solutions)
#     # print(real_solutions)
#     for sol in real_solutions:
#         if current_leaf_length < sol <= L:
#             unique_solution = sol
#             leaf_length_increment = unique_solution - current_leaf_length
#             break
#         elif sol <= current_leaf_length:
#             leaf_length_increment = 0.0
#         elif sol > L:
#             leaf_length_increment = L - current_leaf_length
    
#     return leaf_length_increment



# def search_new_dl(d, l, dS_input, L):
#     found_start = False  
#     first = True
#     key_prev = None
#     value_prev = None
#     for key, value in d.items():
#         if key >= l:
#             found_start = True  
#             if first:
#                 S_start = value
#                 key_prev, value_prev = key, value
#                 first = False
#         if found_start:
#             # print(key, value-S_start)
#             if value-S_start >= dS_input:  
#                 return key_prev, value_prev 
#             else:
#                 key_prev, value_prev = key, value
#     return L, value  # if no match is found --> finished

