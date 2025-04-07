from openalea.mtg import *

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

def add_tiller(g, vid):
    """
    Add a tiller to the MTG.

    :algo:
    - get the rank of vid to define the length of the tiller
    - add the first internode of the tiller in the MTG
    - first connect the different scales (tillers, phytomers, internodes) to the axis
    - then add the sequence of vertices (add_child)
    """
    # Add a new component to the MTG

    scale_id = g.scale(vid)
    axis_id = g.complex(vid)
    rank = g.Rank(vid) + 1  # Number of edge from the root of the axis
    n = len(g.Axis(vid))
    len_tiller  = n - rank  # we remove the parent that do not belong to the tiller 

    #assert(len(g.Axis(vid)) == 10)

    tid = g.add_child(parent=axis_id, edge_type='+', label='Axis')
    vid, tid2 = g.add_child_and_complex(parent=vid, complex=tid, edge_type='+', label='I1')

    for i in range(1, len_tiller):
        vid = g.add_child(parent=vid, edge_type='<', label=f'I{i+1}')

    assert(tid==tid2)
    return g

