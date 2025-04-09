from openalea.mtg import *
from openalea.mtg.traversal import pre_order


def test_axis(n=10):
    """
    Create a MTG as an axis.
    """
    # Create a simple MTG
    g = MTG()

    vid = g.add_component(complex_id=g.root, label='Plant')
    vid = g.add_component(complex_id=vid, label='Axis')
    vid = g.add_component(complex_id=vid, label='I1')
    for i in range(1,n):
        vid = g.add_child(parent=vid, edge_type='<', label=f'I{i+1}')
    
    return g

def add_tiller(g, vid, start_time, phyllochrone=1, tiller_delay=1):
    """
    Add a tiller to the MTG.

    :algo:
    - get the rank of vid to define the length of the tiller
    - add the first internode of the tiller in the MTG
    - first connect the different scales (tillers, phytomers, internodes) to the axis
    - then add the sequence of vertices (add_child)
    """
    # Add a new component to the MTG

    tillers = []
    scale_id = g.scale(vid)
    axis_id = g.complex(vid)
    rank = g.Rank(vid) + 1  # Number of edge from the root of the axis
    n = len(g.Axis(vid))
    len_tiller  = n - rank  # we remove the parent that do not belong to the tiller 

    #assert(len(g.Axis(vid)) == 10)
    prev_time = start_time

    tid = g.add_child(parent=axis_id, edge_type='+', label='Axis')
    vid, tid2 = g.add_child_and_complex(parent=vid, complex=tid, edge_type='+', label='I1')

    # add timing to elements!
    g.node(vid).start_tt = prev_time
    g.node(vid).end_tt = prev_time + phyllochrone
    prev_time += phyllochrone
    tillers.append((vid, prev_time+tiller_delay))


    for i in range(1, len_tiller):
        vid = g.add_child(parent=vid, edge_type='<', label=f'I{i+1}')
        g.node(vid).start_tt = prev_time
        g.node(vid).end_tt = prev_time + phyllochrone
        prev_time += phyllochrone

        tillers.append((vid, prev_time+tiller_delay))

    return tillers

    

def test_full_plant(nb_leaves=10, nb_tillers=10, phyllochrone=1, tiller_delay=1):
    """
    Create a MTG as a full plant.
    """
    # Here we consider that the list is sorted by the time
    tiller_points = []
    # Create a simple MTG
    g = test_axis(n=nb_leaves)
    max_scale = g.max_scale()
    root_id = next(g.component_roots_at_scale_iter(g.root, scale=max_scale))

    prev_time = 0
    for vid in pre_order(g, root_id):
        g.node(vid).start_tt = prev_time
        g.node(vid).end_tt = prev_time+phyllochrone
        prev_time += phyllochrone

        tiller_points.append((vid, prev_time+tiller_delay))

    for i in range(nb_tillers):
        # add a tiller
        vid, time = tiller_points.pop(0)

        new_tillers = add_tiller(g, vid, time, phyllochrone=phyllochrone, tiller_delay=tiller_delay) 

        tiller_points.extend(new_tillers)
        
        tiller_points.sort(key=lambda x: x[1])
        tiller_points = tiller_points[:nb_tillers-i]


    return g
