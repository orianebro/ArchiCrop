from __future__ import annotations

from openalea.mtg.traversal import pre_order2, pre_order2_with_filter
from openalea.plantgl.all import (Color3, Material, Vector3)
from oawidgets.plantgl import *

from .display import build_scene
from .geometry import CerealsContinuousVisitor, CerealsTurtle
from .grow_from_constraint import CerealsVisitorConstrained, compute_nb_of_growing_organs


def thermal_time(g, phyllochron=50.0, plastochron=40.0, leaf_duration=1.6, stem_duration=1.6):
    """
    Add dynamic properties on the mtg to simulate development
    leaf_duration is the phyllochronic time for a leaf to develop from tip appearance to collar appearance
    stem_duration is the phyllochronic time for a stem to develop
    falling_rate (degrees / phyllochron) is the rate at which leaves fall after colar appearance
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
        dtt_leaf = phyllochron * leaf_duration

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

    return g


def mtg_turtle_time(g, time, update_visitor=None):
    """Compute the geometry on each node of the MTG using Turtle geometry.

    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.

    :Example:

        >>> def grow(node, time):

    """

    g.properties()["geometry"] = {}
    g.properties()["_plant_translation"] = {}

    max_scale = g.max_scale()

    cereal_visitor = CerealsContinuousVisitor(False)

    def traverse_with_turtle_time(g, vid, time, visitor=cereal_visitor):
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

        if g.node(vid).start_tt <= time:
            visitor(g, vid, turtle, time)
            # turtle.push()
        # plant_id = g.complex_at_scale(vid, scale=1)

        for v in pre_order2_with_filter(g, vid, None, push_turtle, pop_turtle):
            if v == vid:
                continue
            # Done for the leaves
            if g.node(v).start_tt > time:
                continue
            visitor(g, v, turtle, time)

        # scene = turtle.getScene()
        return g

    for plant_id in g.component_roots_at_scale_iter(g.root, scale=max_scale):
        g = traverse_with_turtle_time(g, plant_id, time)
    return g

#############################################################


def init_visible_length(g):
    
    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))

        for metamer in pre_order2(g, v):
            nm = g.node(metamer)

            if "Stem" in nm.label or "Leaf" in nm.label:
                nm.visble_length = 0

    return g


def mtg_turtle_time_with_constraint(g, time, constraint_plants, update_visitor=None):
    """Compute the geometry on each node of the MTG using Turtle geometry.

    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.

    :Example:

        >>> def grow(node, time):

    """

    g.properties()["geometry"] = {}
    g.properties()["_plant_translation"] = {}

    max_scale = g.max_scale()


    def distribute_constraint_among_organs(g, time, constraint_plants):

        # compute number of growing organs
        nb_of_growing_internodes, nb_of_growing_leaves = compute_nb_of_growing_organs(g, time)
        
        # retrieve constraint for height and LA
        constraint_plants_height = constraint_plants[0]
        constraint_plants_LA = constraint_plants[1]
    
        # set initial total amount of growth to distribute among organs
        constraint_to_distribute_height = constraint_plants_height
        constraint_to_distribute_LA = constraint_plants_LA
        # print(constraint_to_distribute_LA)
    
        # set initial distribution per organ (H: same amount for all organs)
        if nb_of_growing_leaves == 0:
            constraint_on_each_leaf = 0
        else:
            constraint_on_each_leaf = constraint_to_distribute_LA / nb_of_growing_leaves # and internodes = sheath !!!!
            # print("Constraint on each leaf:", constraint_on_each_leaf)

        if nb_of_growing_internodes == 0:
            constraint_on_each_internode = 0
        else:
            constraint_on_each_internode = constraint_to_distribute_height / nb_of_growing_internodes
    
        nb_of_updated_leaves = 0
        nb_of_updated_internodes = 0

        return constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes

    
    constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes = distribute_constraint_among_organs(g, time, constraint_plants)

    cereal_visitor = CerealsVisitorConstrained(False)
    
    
    def traverse_with_turtle_time(g, vid, time, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes, visitor=cereal_visitor):
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

        if g.node(vid).start_tt <= time:
            constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes = visitor(g, vid, turtle, time, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes)
            # turtle.push()
        # plant_id = g.complex_at_scale(vid, scale=1)

        for v in pre_order2_with_filter(g, vid, None, push_turtle, pop_turtle):
            if v == vid:
                continue
            # Done for the leaves
            if g.node(v).start_tt > time:
                continue
            constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes = visitor(g, v, turtle, time, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes)

        # scene = turtle.getScene()
        return g


    for plant_id in g.component_roots_at_scale_iter(g.root, scale=max_scale):
        g = traverse_with_turtle_time(g, plant_id, time, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes)
    return g


def grow_plant(g, time, phyllochron):
    g = thermal_time(g, phyllochron)
    g = mtg_turtle_time(g, time=time)
    return g


def grow_plant_with_constraint(g, time, phyllochron, constraint):
    g = thermal_time(g, phyllochron)
    g = mtg_turtle_time_with_constraint(g, time=time, constraint=constraint)
    return g



""" def display_in_NB(g, time):
    g, scene, nump=display_plant(g, time)
    w=PlantGL(scene, group_by_color=True)
    w.wireframe=True
    return w """

# display_in_NB(g, tt_cum[0])
    
# max_time = max(g.property('end_tt').values())
# interact(display_in_NB, g=fixed(g), time=IntSlider(min=min(tt_cum), max=max(tt_cum), step=100, value=tt_cum[10]))