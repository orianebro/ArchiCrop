from __future__ import annotations

import sympy as sp
from openalea.mtg.traversal import pre_order2

from .geometry import leaf_mesh_for_growth, stem_mesh, CerealsTurtle, CerealsVisitor, addSets, TurtleFrame


def read_columns_from_file(filename, separator=';'):
    '''returns list of putput data through time from columns of STICS output file'''
    with open(filename, 'r') as file:
        columns = []
        for line in file:
            values = line.strip().split(separator) # Strip any extra whitespace and split the line by the separator
            if not columns:
                columns = [[] for _ in values]
            while len(columns) < len(values):
                columns.append([])
            for i, value in enumerate(values):
                columns[i].append(value)
            # Handle lines with fewer values than columns
            for i in range(len(values), len(columns)):
                columns[i].append('')
    return columns


def distribute_constraint_among_plants(constraints_crop, nb_of_plants):
    constraint_crop_height = constraints_crop[0]
    constraint_crop_LAI = constraints_crop[1]
    constraint_plants_height = constraint_crop_height
    constraint_plants_LA = constraint_crop_LAI / nb_of_plants
    constraints_plants = [constraint_plants_height, constraint_plants_LA]
    return constraints_plants


def compute_nb_of_growing_organs(g, time):

    nb_of_growing_internodes = 0
    nb_of_growing_leaves = 0

    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))

        for metamer in pre_order2(g, v):
            
            n = g.node(metamer)
            if n.label.startswith("Leaf") or n.label.startswith("Stem"):
                if n.start_tt <= time < n.end_tt: # and n.visible_length < n.mature_length:  # organ growing   
                    if n.label.startswith("Leaf"):
                        nb_of_growing_leaves += 1 
                    if n.label.startswith("Stem"):
                        nb_of_growing_internodes += 1

    return nb_of_growing_internodes, nb_of_growing_leaves


def compute_leaf_length_increment(constraint_LA, leaf):
    """returns the leaf length of a leaf given a leaf area"""

    # def scaled_leaf_shape(s, L, alpha=-2.3):
    #     beta = -2 * (alpha + np.sqrt(-alpha))
    #     gamma = 2 * np.sqrt(-alpha) + alpha
    #     r = alpha * (s / L) ** 2 + beta * (s / L) + gamma
    #     return r


    L = leaf.mature_length
    current_leaf_length = leaf.visible_length

    # # compute current leaf area
    # if current_leaf_length == 0:
    #     current_blade_area = 0
    # else:
    #     lower_bound = max(L - current_leaf_length, 0.0)
    #     upper_bound = L
    #     current_blade_area, error = quad(
    #         scaled_leaf_shape, lower_bound, upper_bound, args=(L, alpha)
    #     )
    #     current_blade_area = 2 * leaf.shape_max_width * blade_area
    # print(current_blade_area)

    
    def leaf_area_function_of_length(s, L, wl=0.12):
        '''returns the current leaf area (i.e. twice the integral of leaf shape) given the current length s
        
            function obtained this way : 

            # Define the variable and parameters
            s = sp.symbols('s')
            L, wl = sp.symbols('L, wl') 

            # Define the scaled leaf shape function with parameters
            beta = -2 * (-2.3 + math.sqrt(2.3))
            gamma = 2 * math.sqrt(2.3) - 2.3
            r = wl * L * (-2.3 * (s / L) ** 2 + beta * (s / L) + gamma)

            # Reflect the function to the axis x=0.5 
            reflected_r = r.subs(s, -s+L).simplify()

            # Find the indefinite integral (primitive)
            primitive_f = sp.integrate(reflected_r, s)
        
        '''

        return 2 * (1.51657508881031*s**2*wl - 0.766666666666667*s**3*wl/L)
    
    s = sp.symbols('s')
    
    # Define the equation
    equation = constraint_LA + leaf_area_function_of_length(current_leaf_length, L) - leaf_area_function_of_length(s, L)
    
    # Solve the equation, restricting the domain to real numbers
    real_solutions = sp.solveset(equation, s, domain=sp.S.Reals)
    # print(real_solutions)
    for sol in real_solutions:
        if current_leaf_length < sol < L:
            unique_solution = sol
            # break

    leaf_length_increment = unique_solution - current_leaf_length
    
    # print("The unique solution is :", unique_solution)
    # print("The leaf length increment is :", leaf_length_increment)
    # print("The leaf area increment is :", leaf_area_function_of_length(unique_solution, L) - leaf_area_function_of_length(current_leaf_length, L))
    
    return leaf_length_increment


def compute_continuous_element_with_constraint(
    element_node, time, constraint_on_each_leaf, constraint_on_each_internode, 
    constraint_to_distribute_LA, constraint_to_distribute_height, 
    nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes, classic=False
):  # see maybe with *kwarg, **kwds, etc. for time
    """compute geometry of Adel base elements (LeafElement and StemElement)
    element_node should be a mtg node proxy"""
    n = element_node
    geom = None

    # added by Oriane, for continuous growth with constraint
    
    # starting from the lowest (oldest) growing organ
    if n.start_tt <= time < n.end_tt and n.visible_length < n.mature_length:  # organ growing
        # for constrained growth
        
        if n.label.startswith("Leaf"):
            # convert constraint for each organ(t) to leaf length
            leaf_length_increment = compute_leaf_length_increment(constraint_on_each_leaf, n)
            # print(leaf_length_increment, n.mature_length)
            
            if n.visible_length + leaf_length_increment <= n.mature_length:
                n.visible_length += leaf_length_increment
                nb_of_updated_leaves += 1
                constraint_to_distribute_LA -= constraint_on_each_leaf
            else: 
                constraint_to_distribute_LA -= 1.5 * n.wl * (n.mature_length**2 - n.visible_length**2) 
                n.visible_length = n.mature_length
                nb_of_updated_leaves += 1
                constraint_on_each_leaf = constraint_to_distribute_LA / (nb_of_growing_leaves - nb_of_updated_leaves)
            n.grow = True
            
        elif n.label.startswith("Stem"):
            ''' 
            h = n.visible_length
            radius = n.diameter / 2
            sheath_area = 2 * np.pi * radius * h
            '''
            if n.visible_length + constraint_on_each_internode <= n.mature_length:
                n.visible_length += constraint_on_each_internode
                nb_of_updated_internodes += 1
                constraint_to_distribute_height -= constraint_on_each_internode
            else: 
                constraint_to_distribute_height -= (n.mature_length - n.visible_length)
                n.visible_length = n.mature_length
                nb_of_updated_internodes += 1
                constraint_on_each_internode = constraint_to_distribute_height / (nb_of_growing_internodes - nb_of_updated_internodes)
            n.grow = True

    elif time < n.start_tt:  # organ not yet appeared
        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            n.visible_length = 0.0
            n.grow = False

    elif time >= n.end_tt and n.visible_length < n.mature_length:  # organ reaches maturity below potential 
        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            n.visible_length = n.visible_length
            n.grow = True

    elif time >= n.end_tt:  # organ reaches maturity at potential
        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            n.visible_length = n.mature_length
            n.grow = True

    
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
                    stem_diameter=n.stem_diameter,
                )
            if n.lrolled > 0:
                rolled = stem_mesh(
                    n.lrolled, n.lrolled, n.d_rolled, n.d_rolled, classic
                )
                if geom is None:
                    geom = rolled
                else:
                    geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
    elif n.label.startswith("Stem"):  # stem element
        geom = stem_mesh(n.length, n.visible_length, n.diameter, n.diameter, classic)

    return geom, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes
        
    

class CerealsVisitorConstrained(CerealsVisitor):
    def __init__(self, classic):
        super().__init__(classic)

    def __call__(self, g, v, turtle, time, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes):
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
            mesh, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes = compute_continuous_element_with_constraint(n, time, constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes, self.classic)
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

        return constraint_on_each_leaf, constraint_on_each_internode, constraint_to_distribute_LA, constraint_to_distribute_height, nb_of_growing_leaves, nb_of_updated_leaves, nb_of_growing_internodes, nb_of_updated_internodes


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