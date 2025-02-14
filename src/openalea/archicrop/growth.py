from __future__ import annotations

import math
import numpy as np
import sympy as sp
from scipy.interpolate import splev
from openalea.mtg.traversal import pre_order2, pre_order2_with_filter

from .geometry import CerealsTurtle, CerealsVisitor, addSets, TurtleFrame, leaf_mesh_for_growth, stem_mesh
from .plant_shape import shape_to_surface




def thermal_time(g, phyllochron, plastochron, leaf_duration, stem_duration): # durations 1.6
    """
    Add dynamic properties on the mtg to simulate development

    :param g: MTG, MTG of a plant
    :param phyllochron: float, phyllochron, i.e. internode appearance rate (in °C.day/internode)
    :param plastochron: float, plastochron, i.e. leaf appearance rate (in °C.day/leaf)
    :param leaf_duration: float, phyllochronic time for a leaf to develop from tip appearance to collar appearance (/phyllochron)
    :param stem_duration: float, phyllochronic time for a stem to develop from base to top (/phyllochron)
    # falling_rate (degrees / phyllochron) is the rate at which leaves fall after colar appearance

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

                nm.senescence = nm.end_tt + 1000 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return g


def init_params_for_growth(g):
    """Initialize leaf and stem lengths and stem diameters"""

    for k in g.properties()["visible_length"].keys():
        g.properties()["visible_length"][k] = 0.0

    
    # axes = g.vertices(scale=1)
    # metamer_scale = g.max_scale()

    # for axis in axes:
    #     v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))

    #     for metamer in pre_order2(g, v):
    #         nm = g.node(metamer)

    #         if nm.label.startswith("Leaf") or nm.label.startswith("Stem"):
    #             nm.visible_length = 0.0
                # nm.stem_diameter = 0.0
                
    # for id in g.vertices():
    #     print(g[id])
    # print(" ")

    return g


def init_growth_dict():
    """Initialize the dictionnary of variables recording growth process"""

    return {"LA_for_each_leaf": 0.0,
            "height_for_each_internode": 0.0,
            "LA_to_distribute": 0.0,
            "height_to_distribute": 0.0,
            "nb_of_growing_leaves": 0,
            "nb_of_growing_internodes": 0,
            "nb_of_updated_leaves": 0,
            "nb_of_updated_internodes": 0}


def compute_nb_of_growing_organs(g, time):
    """Compute the number of growing organs for a given time step"""

    pot_length_growing_internodes = []
    pot_area_growing_leaves = []

    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))

        for metamer in pre_order2(g, v):
            
            n = g.node(metamer)
            if n.label.startswith("Leaf"):
                if n.start_tt <= time < n.end_tt: # and n.visible_length < n.mature_length:  # organ growing   
                    pot_area_growing_leaves.append(n.leaf_area)

            if n.label.startswith("Stem"):
                if n.start_tt <= time < n.end_tt: # and n.visible_length < n.mature_length:  # organ growing  
                    pot_length_growing_internodes.append(n.mature_length) 

    return pot_length_growing_internodes, pot_area_growing_leaves


def distribute_constraint_among_organs(g, time, increments, growth):
    """Fill the dictionnary of variables recording growth process for a give time step"""

    # compute number of growing organs
    pot_length_growing_internodes, pot_area_growing_leaves = compute_nb_of_growing_organs(g, time)
    nb_of_growing_internodes, nb_of_growing_leaves = len(pot_length_growing_internodes), len(pot_area_growing_leaves)

    # set initial total amount of growth to distribute among organs
    # height_to_distribute = increments["Height increment"] + growth["height_to_distribute"]
    # LA_to_distribute = increments["Leaf area increment"] + growth["LA_to_distribute"]
    height_to_distribute = increments["Height increment"]
    LA_to_distribute = increments["Leaf area increment"]

    def equal_distribution(ds, ng):
        '''return initial distribution per organ (H: same amount for all organs)'''
        return [ds/g] * len(ng)
    
    def demand_related_distribution(increment, pot_growing_organs):
        '''return distribution of increment per organ proportionnal to potential length ou area'''
        sum_growing_organs = sum(pot_growing_organs)
        return [increment*p/sum_growing_organs for p in pot_growing_organs]
    
    if nb_of_growing_leaves != 0:
        # LA_for_each_leaf = equal_distribution(LA_to_distribute, nb_of_growing_leaves)
        LA_for_each_leaf = demand_related_distribution(LA_to_distribute, pot_area_growing_leaves)
    else:
        LA_for_each_leaf = 0.0

    if nb_of_growing_internodes != 0:
        # LA_for_each_leaf = equal_distribution(height_to_distribute, nb_of_growing_internodes)
        height_for_each_internode = demand_related_distribution(height_to_distribute, pot_length_growing_internodes)
    else:
        height_for_each_internode = 0.0

    nb_of_updated_leaves = 0
    nb_of_updated_internodes = 0

    growth = {"LA_for_each_leaf": LA_for_each_leaf,
              "height_for_each_internode": height_for_each_internode,
              "LA_to_distribute": LA_to_distribute,
              "height_to_distribute": height_to_distribute,
              "nb_of_growing_leaves": nb_of_growing_leaves,
              "nb_of_growing_internodes": nb_of_growing_internodes,
              "nb_of_updated_leaves": nb_of_updated_leaves,
              "nb_of_updated_internodes": nb_of_updated_internodes}

    return growth


def leaf_area_function_of_length(l, L, wl, alpha=-2.3):
    '''returns the current leaf area (i.e. twice the integral of leaf shape) given the current length s
    
        function obtained this way : 

        # Define the variable and parameters
        s = sympy.symbols('l')
        L, wl, alpha = sympy.symbols('L, wl, alpha') 

        # Define the scaled leaf shape function with parameters
        # alpha = -2.3
        beta = -2 * (alpha + sqrt(-alpha))
        gamma = 2 * sqrt(-alpha) + alpha
        r = wl * L (alpha * s**2 + beta * s + gamma)

        # Reflect the function to the axis x=0.5 
        reflected_r = r.subs(s, -s+L).simplify()

        # Find the indefinite integral (primitive)
        primitive_f = sympy.integrate(reflected_r, s)

        print(2*primitive_f)
    
    '''

    return 2*l**2*wl*math.sqrt(-alpha) + 2*alpha*l**3*wl/(3*L)


def compute_leaf_length_increment(LA_for_each_leaf, leaf):
    """returns the leaf length given a leaf area"""

    L = leaf.mature_length
    current_leaf_length = leaf.visible_length
    wl = leaf.shape_max_width / leaf.mature_length
    
    s = sp.symbols('s')
    
    # Define the equation
    equation = LA_for_each_leaf + leaf_area_function_of_length(current_leaf_length, L, wl) - leaf_area_function_of_length(s, L, wl)
    # print ("LA_for_each_leaf", LA_for_each_leaf)

    # Solve the equation, restricting the domain to real numbers
    real_solutions = sp.solveset(equation, s, domain=sp.S.Reals)
    # print("Sol", real_solutions)
    # print(real_solutions)
    for sol in real_solutions:
        if current_leaf_length < sol <= L:
            unique_solution = sol
            leaf_length_increment = unique_solution - current_leaf_length
            break
        elif sol <= current_leaf_length:
            leaf_length_increment = 0.0
        elif sol > L:
            leaf_length_increment = L - current_leaf_length
    
    return leaf_length_increment



def search_new_dl(d, l, dS_input, L):
    found_start = False  
    first = True
    key_prev = None
    value_prev = None
    for key, value in d.items():
        if key >= l:
            found_start = True  
            if first:
                S_start = value
                key_prev, value_prev = key, value
                first = False
        if found_start:
            # print(key, value-S_start)
            if value-S_start >= dS_input:  
                return key_prev, value_prev 
            else:
                key_prev, value_prev = key, value
    return L, value  # if no match is found --> finished



def compute_continuous_element_with_constraint(
    element_node, time, growth, classic=False
):  # see maybe with *kwarg, **kwds, etc. for time
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""
    n = element_node
    geom = None

    # added by Oriane, for continuous constrained growth 

    LA_for_each_leaf = growth["LA_for_each_leaf"]
    height_for_each_internode = growth["height_for_each_internode"]
    LA_to_distribute = growth["LA_to_distribute"]
    height_to_distribute = growth["height_to_distribute"]
    nb_of_growing_leaves = growth["nb_of_growing_leaves"]
    nb_of_growing_internodes = growth["nb_of_growing_internodes"]
    nb_of_updated_leaves = growth["nb_of_updated_leaves"]
    nb_of_updated_internodes = growth["nb_of_updated_internodes"]


    if n.label.startswith("Leaf") or n.label.startswith("Stem"):

        if time < n.start_tt:  # organ not yet appeared
            n.visible_length = 0.0
            # n.stem_diameter = 0.0
            n.grow = False
            n.lengths.append(0.0)
    
        elif n.start_tt <= time:
            n.grow = True
            # update age of organ
            n.age = time - n.start_tt
        
    
            # starting from the lowest (oldest) growing organ
            if time <= n.end_tt: # if n.start_tt <= time < n.end_tt: # organ growing

                # n.stem_diameter = 1.2

                if n.visible_length < n.mature_length:  

                    if n.label.startswith("Leaf"):
                        # convert constraint on area for each leaf(t) to leaf length
                        LA_for_this_leaf = LA_for_each_leaf[nb_of_updated_leaves]
                        
                        if n.visible_leaf_area + LA_for_this_leaf <= n.leaf_area:
                            n.visible_leaf_area += LA_for_this_leaf
                            # leaf_length_increment = compute_leaf_length_increment(LA_for_this_leaf, n)
                            # n.visible_length += leaf_length_increment
                            percent_of_progression_in_area = n.visible_leaf_area / n.leaf_area
                            tck = shape_to_surface(n.shape, n.wl) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            n.visible_length = float(splev(x=percent_of_progression_in_area, tck=tck)) * n.mature_length
                        else: 
                            LA_rest = LA_for_this_leaf - (n.leaf_area - n.visible_leaf_area)
                            n.visible_leaf_area = n.leaf_area
                            next_leaves_to_update = len(LA_for_each_leaf[nb_of_updated_leaves+1:])
                            for i in range(next_leaves_to_update):
                                LA_for_each_leaf[nb_of_updated_leaves+1+i] += LA_rest/next_leaves_to_update
                            n.visible_length = n.mature_length
                        nb_of_updated_leaves += 1
                        n.lengths.append(n.visible_length)
                            # if nb_of_growing_leaves > nb_of_updated_leaves:
                            #     LA_for_each_leaf = LA_to_distribute / (nb_of_growing_leaves - nb_of_updated_leaves)
                            # else:
                            #     LA_for_each_leaf = 0.0

                        
                
                    elif n.label.startswith("Stem"):

                        height_for_this_internode = height_for_each_internode[nb_of_updated_internodes]

                        if n.visible_length + height_for_this_internode <= n.mature_length:
                            n.visible_length += height_for_this_internode
                        else: 
                            height_rest = height_for_this_internode - (n.mature_length - n.visible_length)
                            next_internodes_to_update = len(height_for_each_internode[nb_of_updated_internodes+1:])
                            for i in range(next_internodes_to_update):
                                height_for_each_internode[nb_of_updated_internodes+1+i] += height_rest/next_internodes_to_update
                            n.visible_length = n.mature_length
                        nb_of_updated_internodes += 1
                        n.lengths.append(n.visible_length)
                            # if nb_of_growing_internodes > nb_of_updated_internodes:
                            #     height_for_each_internode = height_to_distribute / (nb_of_growing_internodes - nb_of_updated_internodes)

                # else:
                #     n.visible_length = n.mature_length

                
            elif time > n.end_tt:
                # n.stem_diameter = min(1.2 + (n.mature_stem_diameter - 1.2)*0.001*time, n.mature_stem_diameter)
                if n.visible_length < n.mature_length:  # organ reaches maturity below potential 
                    n.visible_length = n.visible_length
                else:  # organ reaches maturity at potential
                    n.visible_length = n.mature_length
                n.lengths.append(n.visible_length)

                #
                #
                #
                #

                if n.label.startswith("Leaf"):
                    if n.senescence <= time:
                        # n.visible_length = 0.0
                        n.grow = False
                
                #
                #
                #
                #

    # update growth dict
    growth["LA_for_each_leaf"] = LA_for_each_leaf
    growth["height_for_each_internode"] = height_for_each_internode
    growth["LA_to_distribute"] = LA_to_distribute
    growth["height_to_distribute"] = height_to_distribute
    growth["nb_of_growing_leaves"] = nb_of_growing_leaves
    growth["nb_of_growing_internodes"] = nb_of_growing_internodes
    growth["nb_of_updated_leaves"] = nb_of_updated_leaves
    growth["nb_of_updated_internodes"] = nb_of_updated_internodes

    
    
    if n.label.startswith("Leaf"):  # leaf element
        if n.visible_length > 0.0001:  # filter less than 0.001 mm leaves
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
        

class CerealsVisitorConstrained(CerealsVisitor):
    def __init__(self, classic):
        super().__init__(classic)

    def __call__(self, g, v, turtle, time, growth):
        # 1. retrieve the node
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
            # nb_of_growing_internodes, nb_of_growing_leaves = compute_nb_of_growing_organs(g, time)
            mesh, growth = compute_continuous_element_with_constraint(n, time, growth, self.classic)
            if mesh:  # To DO : reset to None if calculated so ?
                n.geometry = turtle.transform(mesh)
                n.anchor_point = turtle.getPosition()

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


def mtg_turtle_time_with_constraint(g, time, increments, growth, update_visitor=None):
    """Compute the geometry on each node of the MTG using Turtle geometry.

    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.

    :Example:

        >>> def grow(node, time):

    """

    g.properties()["geometry"] = {}
    g.properties()["_plant_translation"] = {}

    max_scale = g.max_scale()

    growth = distribute_constraint_among_organs(g, time, increments, growth)

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
    return g, growth


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



