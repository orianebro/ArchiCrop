from __future__ import annotations

from openalea.mtg.traversal import pre_order2, pre_order2_with_filter
from openalea.plantgl.all import (
    Color3,
    Material,
)

from .display import build_scene
from .geometry import CerealsContinuousVisitor, CerealsTurtle


def thermal_time(g, phyllochron=110., leaf_duration=1.6, stem_duration=1.6):
    """
    Add dynamic properties on the mtg to simulate developpement
    leaf_duration is the phyllochronic time for a leaf to develop from tip appearance to collar appearance
    stem_duration is the phyllochronic time for a stem to develop
    falling_rate (degrees / phyllochron) is the rate at which leaves fall after colar appearance
    """

    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        tt = 0
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))
        # stem_ids = g.Trunk(v)
        # nb_stems = len(stem_ids)
        # nb_sectors = 1
        dtt_stem = phyllochron*stem_duration
        dtt_leaf = phyllochron*leaf_duration
        
        for metamer in pre_order2(g, v):
            nm = g.node(metamer)

            if 'Stem' in nm.label:
                nm.start_tt = tt
                nm.end_tt = tt + dtt_stem
            elif 'Leaf' in nm.label:
                nm.start_tt = tt
                nm.end_tt = tt + dtt_leaf
                tt += phyllochron

    return g



def mtg_turtle_time(g, time, update_visitor=None):
    ''' Compute the geometry on each node of the MTG using Turtle geometry. 
    
    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.
    
    :Example:

        >>> def grow(node, time):
                
    '''

    g.properties()['geometry'] = {}
    g.properties()['_plant_translation'] = {}

    max_scale = g.max_scale()

    cereal_visitor = CerealsContinuousVisitor(False)

    def traverse_with_turtle_time(g, vid, time, visitor=cereal_visitor):
        turtle = CerealsTurtle()
        def push_turtle(v):
            n = g.node(v)
            #if 'Leaf' in n.label:
                #    return False
            try:
                start_tt = n.start_tt
                if start_tt > time:
                    return False
            except: 
                pass
            if g.edge_type(v) == '+':
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
            if g.edge_type(v) == '+':  # noqa: RET503
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



def grow_plant(g, time):
    g = thermal_time(g)
    g = mtg_turtle_time(g, time=time)
    return g


def grow_plant_and_display(g, time):
    g = grow_plant(g=g, time=time)
    # Build and display scene
    nice_green=Color3((50,100,0))
    scene, nump = build_scene(g, 
                              leaf_material=Material(nice_green), 
                              stem_material=Material(nice_green))
    return g, scene, nump