from __future__ import annotations

from openalea.mtg import *  # noqa: F403
from openalea.mtg.traversal import pre_order
from openalea.mtg.turtle import traverse_with_turtle


def main_axis(n=10):
    """
    Create a MTG as an axis.
    """
    # Create a simple MTG
    g = MTG()  # noqa: F405

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
    scale_id = g.scale(vid)  # noqa: F841
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

    

def generate_full_plant(nb_leaves=10, nb_tillers=10, phyllochrone=1, tiller_delay=1):
    """
    Create a MTG as a full plant.
    """
    # Here we consider that the list is sorted by the time
    tiller_points = []
    # Create a simple MTG
    g = main_axis(n=nb_leaves)
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

# tests
# g=generate_full_plant(nb_leaves=10, nb_tillers=1)
# g=generate_full_plant(nb_leaves=10, nb_tillers=2)
# g=generate_full_plant(nb_leaves=15, nb_tillers=40)
# g.display()

def nb_ramifs(g):
    max_scale = g.max_scale()
    ramifs = [v for v in g.vertices_iter(scale=max_scale-1) if g.edge_type(v)=='+']
    return len(ramifs)



def test_main_axis():
    """
    Test the main_axis function.
    """
    g = generate_full_plant(nb_leaves=10, nb_tillers=0)
    assert len(g) == 13
    assert g.root == 0
    assert g.scale(3) == 3

    assert nb_ramifs(g) == 0

def test_t1():
    g=generate_full_plant(nb_leaves=10, nb_tillers=1)
    assert nb_ramifs(g) == 1

def test_t5():
    g=generate_full_plant(nb_leaves=10, nb_tillers=5)
    assert nb_ramifs(g) == 5

def test_t40():
    g=generate_full_plant(nb_leaves=20, nb_tillers=40)
    assert nb_ramifs(g) == 40


def visitor(g, v, turtle):
    if g.edge_type(v) == '+':
        angle = 60 if g.order(v) == 1 else 30
        turtle.down(angle)
    turtle.setId(v)
    turtle.F(1)
    turtle.rollL()

def plot(g):
    # Viewer.display(plot(g))
    scene = traverse_with_turtle(g, 3, visitor=visitor, gc=True)
    return scene  # noqa: RET504


'''
# Print leaf area of all leaves, axis by axis
la_per_axis = {}
i = 0
for vid in g.vertices(scale=2):  # scale=2 corresponds to axes
    leaf_areas = []
    axis_label = g.parent(vid)
    print(f"Axis {axis_label}: {len(g.components(vid))/2} leaves")
    for leaf_vid in g.components(vid):
        if g.label(leaf_vid) == "Leaf":
            la = g.property("leaf_area")[leaf_vid]
            # la = g.property("mature_length")[leaf_vid]
            print(f"  Leaf {leaf_vid}: area = {la}")
            leaf_areas.append(la)
    la_per_axis[vid] = leaf_areas
    i += 1

import matplotlib.pyplot as plt
total_area = 0.0
for key, value in la_per_axis.items():
    total_area += sum(value)
    plt.plot(range(1,len(value)+1), value, label=f"{key}")
print(total_area)
print(leaf_area)
plt.legend()
plt.show()

axes = g.vertices(scale=2)
axis_orders = {axis: g.order(axis) for axis in axes}
for vid in g.vertices(scale=2):  # scale=2 corresponds to axes
    axis_label = g.label(vid)
    print(f"Axis {axis_label}:")
    for leaf_vid in g.components(vid):
        if g.label(leaf_vid) == "Leaf":
            leaf_area = g.property("leaf_area")[leaf_vid]
            print(f"  Leaf {leaf_vid}: area = {leaf_area}")

import oawidgets.mtg
oawidgets.mtg.plot(g)
'''