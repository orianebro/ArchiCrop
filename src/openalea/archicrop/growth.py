from __future__ import annotations

from scipy.interpolate import splev

from openalea.mtg.traversal import pre_order2_with_filter

from .cereal_plant import CerealsVisitor
from .geometry import compute_growing_organ_geometry
from .turtle import CerealsTurtle


def init_visible_variables(g, daily_dynamics):
    """Initialize leaf and stem visible lengths"""

    for k in g.properties()["visible_length"]:
        g.properties()["visible_length"][k] = 0.0
    for k in g.properties()["leaf_lengths"]:
        g.properties()["leaf_lengths"][k] = [0.0]*len(daily_dynamics)
        g.properties()["leaf_areas"][k] = [0.0]*len(daily_dynamics)
        g.properties()["senescent_lengths"][k] = [0.0]*len(daily_dynamics)
    # for k in g.properties()["visible_leaf_area"]:
    #     g.properties()["visible_leaf_area"][k] = 0.0
    for k in g.properties()["leaf_lengths"]:
        g.properties()["stem_lengths"][k] = [0.0]*len(daily_dynamics)
    for k in g.properties()["stem_diameter"]:
        g.properties()["stem_diameter"][k] = 0.0

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

    for vid, la in g.property("leaf_area").items(): 
        n = g.node(vid)
        if n.start_tt <= time <= n.end_tt or prev_time < n.end_tt < time:
            growing_leaves[vid] = {"potential": la, "visible": n.visible_leaf_area}
        if n.senescence <= time and n.visible_leaf_area > n.senescent_area: # not n.dead: # and n.srt > 0: 
            senescing_leaves[vid] = {"potential": n.visible_leaf_area, "visible": n.senescent_area}
        # print(n.senescence <= time)
        # print(n.visible_leaf_area > n.senescent_area)
        # print(senescing_leaves)

    for vid, ml in g.property("mature_length").items(): 
        n = g.node(vid)
        if n.label.startswith("Stem") and (n.start_tt <= time <= n.end_tt or prev_time < n.end_tt < time):
            growing_internodes[vid] = {"potential": ml, "visible": n.visible_length}

    return growing_internodes, growing_leaves, senescing_leaves

'''
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
'''

# rename to grow_organs
def distribute_to_potential(g, growing_organs, day, time, time_increment, increment_to_distribute, distribution_function, grow_function, to_print=False):
    """Distribute increment among growing organs up to potential of each organ"""

    increment_for_each_organ = dict.fromkeys(growing_organs, 0.0)

    
    # if g and g.node(list(growing_organs.keys())[0]).label.startswith("Leaf"):
    #     print(f'To distribute {increment_to_distribute}')
    #     visible_leaf_area = g.property("visible_leaf_area")
    #     print(visible_leaf_area)
    #     print(growing_organs)
    #     surface1= sum(visible_leaf_area.values())
    #     print(f'Plant Leaf Area {surface1}')
        # surface_growing1 = sum(visible_leaf_area[vid] for vid in growing_organs)
        # to_distribute = increment_to_distribute

    while len(growing_organs) > 0 and increment_to_distribute > 1e-5: 
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
        
    for vid, incr in increment_for_each_organ.items():
        n = g.node(vid)
        n.grow = True
        grow_function(n, incr, day, time_increment)

        n_stem = n.parent() if n.label.startswith("Leaf") else n
        n.stem_diameter = min(n.mature_stem_diameter/2 * (1+0.5*(time - n_stem.start_tt) / (n_stem.end_tt - n_stem.start_tt)), n.mature_stem_diameter) 
        
        if n.label.startswith("Leaf"):
            n.inclination = min(1.5*(time - n.start_tt) / (n.end_tt - n.start_tt), 1)

    # if g and g.node(list(growing_organs.keys())[0]).label.startswith("Leaf"):
    #     surface_growing2 = sum(growing_organs[vid]["visible"] for vid in growing_organs)
    #     delta = surface_growing2 - surface_growing1 
    #     if abs(delta-to_distribute)>1e-3 or to_print:
    #         print(f'FAILURE : {delta} and {to_distribute}')
    #         print(f'real delta {delta - sum(increment_for_each_organ.values())}')
        
    #print(f'Plant Leaf Area {sum(g.properties("visible_leaf_area"))}')

    # return increment_for_each_organ



def distribute_among_organs(g, day, daily_dynamics, rate=True):
    """Fill the dictionnary of variables recording growth process at a given time"""

    time = daily_dynamics["Thermal time"]
    thermal_time_incr = daily_dynamics["Thermal time increment"]

    # compute number of growing organs
    growing_internodes, growing_leaves, senescing_leaves = get_growing_and_senescing_organs_potential_visible(g, time, thermal_time_incr)

    # compute number of senscent organs
    # senescing_leaves = get_senescing_organs(g, time)

    # set amount of growth to distribute among organs

    # to_print = 551 <time<552
    if daily_dynamics["Leaf area increment"] > 0.0 and len(growing_leaves) > 0:
        update_cereal_leaf_growth = update_cereal_leaf_growth_rate if rate else update_cereal_leaf_growth_area
        distribute_to_potential(g=g, 
                                growing_organs=growing_leaves, 
                                day=day,
                                time=time,
                                time_increment=thermal_time_incr,
                                increment_to_distribute=daily_dynamics["Leaf area increment"], 
                                distribution_function=demand_dist, 
                                grow_function=update_cereal_leaf_growth)
                                # to_print=to_print)
    # else:
    #     LA_for_each_leaf = dict.fromkeys(growing_leaves, 0.0)
    
    # axes = {}
    # for vid in growing_internodes:
    #     axes.setdefault(g.complex(vid), []).append(vid)

    if daily_dynamics["Height increment"] > 0.0 and len(growing_internodes) > 0:
        update_stem_growth = update_stem_growth_rate if rate else update_stem_growth_height
        distribute_to_potential(g=g, 
                                growing_organs=growing_internodes, 
                                day=day,
                                time=time,
                                time_increment=thermal_time_incr,
                                increment_to_distribute=daily_dynamics["Height increment"], 
                                distribution_function=demand_dist,
                                grow_function=update_stem_growth)
        # height_for_each_internode = {}
        # for aid in axes:  
        #     growing_inter = {vid:growing_internodes[vid]  for vid in axes[aid]}
        #     height_for_each_axis = distribute_to_potential(g=g, 
        #                                                    growing_organs=growing_inter, 
        #                                                    day=day,
        #                                                    time=time,
        #                                                    time_increment=thermal_time_incr,
        #                                                    increment_to_distribute=daily_dynamics["Height increment"], 
        #                                                    distribution_function=demand_dist,
        #                                                    grow_function=update_stem_growth)
        #     height_for_each_internode.update(height_for_each_axis)
    # else:
    #     height_for_each_internode = dict.fromkeys(growing_internodes, 0.0)

    if daily_dynamics["Senescent leaf area increment"] > 0.0: # and len(senescing_leaves) > 0:
        distribute_to_potential(g=g, 
                                growing_organs=senescing_leaves, 
                                day=day,
                                time=time,
                                time_increment=thermal_time_incr,
                                increment_to_distribute=daily_dynamics["Senescent leaf area increment"], 
                                distribution_function=equal_dist,
                                grow_function=update_leaf_senescence_area)

    # else:
    #     sen_LA_for_each_leaf = dict.fromkeys(senescing_leaves, 0.0)

    # return {"thermal_time_increment": thermal_time_incr,
    #         # "height_to_distribute": height_to_distribute,
    #         "height_for_each_internode": height_for_each_internode,
    #         # "LA_to_distribute": LA_to_distribute,
    #         "LA_for_each_leaf": LA_for_each_leaf,
    #         # "sen_LA_to_distribute": sen_LA_to_distribute,
    #         "sen_LA_for_each_leaf": sen_LA_for_each_leaf,
    #         "growing_internodes": growing_internodes,
    #         "growing_leaves": growing_leaves,
    #         "senescing_leaves": senescing_leaves}

'''
def compute_organ_growth(vid, element_node, time, growth, rate=True):
    """Compute the gain or loss in length or area of an organ at a given time."""
    
    n = element_node

    # # is there LA to update for this node ?
    # if vid in growth["LA_for_each_leaf"]:
    #     LA_for_this_leaf = growth["LA_for_each_leaf"][vid]
    #     if rate:
    #         update_cereal_leaf_growth_rate(n, LA_for_this_leaf, growth["thermal_time_increment"])
    #     else:
    #         update_cereal_leaf_growth_area(n, LA_for_this_leaf)


    if n.label.startswith("Leaf") or n.label.startswith("Stem"):
    
        if time < n.start_tt:
            n.visible_length = 0.0
            n.grow = False
            if n.label.startswith("Leaf"):
                n.leaf_lengths.append(0.0)
                n.senescent_lengths.append(0.0)
            elif n.label.startswith("Stem"):
                n.stem_lengths.append(0.0)
        if n.start_tt <= time:
            n.age = time - n.start_tt
        
            if time <= n.end_tt: 
                n.grow = True
                if n.visible_length < n.mature_length:
                    if n.label.startswith("Leaf") and growth["LA_to_distribute"] > 0.0:
                        LA_for_this_leaf = growth["LA_for_each_leaf"][vid]
                        if rate:
                            update_cereal_leaf_growth_rate(n, LA_for_this_leaf, growth["thermal_time_increment"])
                        else:
                            update_cereal_leaf_growth_area(n, LA_for_this_leaf)
                    if n.label.startswith("Stem") and growth["height_to_distribute"] > 0.0:
                        height_for_this_internode = growth["height_for_each_internode"][vid]
                        if rate:
                            update_stem_growth_rate(n, height_for_this_internode, growth["thermal_time_increment"]) 
                        else:
                            update_stem_growth_height(n, height_for_this_internode)

            elif time > n.end_tt:
                if n.visible_length < n.mature_length:
                    if n.label.startswith("Leaf") and growth["LA_to_distribute"] > 0.0 and vid in growth["LA_for_each_leaf"]:
                        LA_for_this_leaf = growth["LA_for_each_leaf"][vid]
                        if rate:
                            update_cereal_leaf_growth_rate(n, LA_for_this_leaf, growth["thermal_time_increment"])
                        else:
                            update_cereal_leaf_growth_area(n, LA_for_this_leaf)
                    if n.label.startswith("Stem") and growth["height_to_distribute"] > 0.0 and vid in growth["height_for_each_internode"]:
                        height_for_this_internode = growth["height_for_each_internode"][vid]
                        if rate:
                            update_stem_growth_rate(n, height_for_this_internode, growth["thermal_time_increment"]) 
                        else:
                            update_stem_growth_height(n, height_for_this_internode)

            if time > n.end_tt:
                if n.label.startswith("Leaf") and growth["sen_LA_to_distribute"] > 0.0:
                    if n.senescence <= time and not n.dead and vid in growth["sen_LA_for_each_leaf"]:
                        sen_LA_for_this_leaf = growth["sen_LA_for_each_leaf"][vid]
                        update_leaf_senescence_area(n, sen_LA_for_this_leaf)
                    else:
                        n.senescent_lengths.append(n.visible_length)
                elif n.label.startswith("Leaf"):
                    n.leaf_lengths.append(n.visible_length)
                    n.senescent_lengths.append(n.senescent_length)
                elif n.label.startswith("Stem"):
                    n.stem_lengths.append(n.visible_length)

        # Evolution of organ geometry with age
        # TODO: take tt from stem of same rank for a leaf

        n_stem = n.parent() if n.label.startswith("Leaf") else n
        n.stem_diameter = min(n.mature_stem_diameter/2 * (1+0.5*(time - n_stem.start_tt) / (n_stem.end_tt - n_stem.start_tt)), n.mature_stem_diameter) 
        
        if n.label.startswith("Leaf"):
            n.inclination = min(1.5*(time - n.start_tt) / (n.end_tt - n.start_tt), 1)


    return growth
'''

def update_cereal_leaf_growth_area(n, LA_for_this_leaf, day, time_increment=None):
    """Update the visible leaf area and length of a cereal leaf, 
    under a potential leaf area."""
    if n.visible_leaf_area + LA_for_this_leaf < n.leaf_area:
        n.visible_leaf_area += LA_for_this_leaf
        relative_visible_area = n.visible_leaf_area / n.leaf_area
        n.visible_length = float(splev(x=relative_visible_area, tck=n.tck)) * n.mature_length
    else:
        n.visible_leaf_area = n.leaf_area
        n.visible_length = n.mature_length
    # print("visible length", n.visible_length)

    for d in range(day-1, len(n.leaf_areas)):
        n.leaf_areas[d] = n.visible_leaf_area
        n.leaf_lengths[d] = n.visible_length
    

def update_cereal_leaf_growth_rate(n, LA_for_this_leaf, day, time_increment):
    """Update the visible leaf area and length of a cereal leaf, 
    under a potential growth rate."""
    rate = LA_for_this_leaf / time_increment
    if rate <= n.potential_growth_rate:
        n.visible_leaf_area += LA_for_this_leaf
    else:
        new_LA_for_this_leaf = n.potential_growth_rate * time_increment
        n.visible_leaf_area += new_LA_for_this_leaf
    relative_visible_area = n.visible_leaf_area / n.leaf_area
    n.visible_length = float(splev(x=relative_visible_area, tck=n.tck)) * n.mature_length

    for d in range(day-1, len(n.leaf_areas)):
        n.leaf_areas[d] = n.visible_leaf_area
        n.leaf_lengths[d] = n.visible_length


def update_stem_growth_height(n, height_for_this_internode, day, time_increment=None):
    """Update the visible length of a stem internode, 
    under a potential height."""
    if n.visible_length + height_for_this_internode <= n.mature_length:
        n.visible_length += height_for_this_internode
    else:
        n.visible_length = n.mature_length

    for d in range(day-1, len(n.stem_lengths)):
        n.stem_lengths[d] = n.visible_length

def update_stem_growth_rate(n, height_for_this_internode, day, time_increment):
    """Update the visible length of a stem internode, 
    under a potential growth rate."""
    rate = height_for_this_internode / time_increment
    if rate <= n.potential_growth_rate:
        n.visible_length += height_for_this_internode
    else:
        new_height_for_this_internode = n.potential_growth_rate * time_increment
        n.visible_length += new_height_for_this_internode

    for d in range(day-1, len(n.stem_lengths)):
        n.stem_lengths[d] = n.visible_length


def update_leaf_senescence_area(n, sen_LA_for_this_leaf, day, time_increment=None):
    """Update the senescent area and length of a cereal leaf, 
    that cannot be greater than the visible leaf area."""
    n.grow = False
    if n.senescent_area + sen_LA_for_this_leaf < n.visible_leaf_area:
        n.senescent_area += sen_LA_for_this_leaf
        relative_senescent_area = n.senescent_area / n.visible_leaf_area
        n.srt = 1 - float(splev(x=relative_senescent_area, tck=n.tck))
        n.senescent_length = (1 - n.srt) * n.visible_length
    else:
        n.senescent_area = n.visible_leaf_area
        n.senescent_length = n.visible_length
        n.dead = True

    # print("senescent length", n.senescent_length)

    for d in range(day, len(n.senescent_lengths)):
        n.senescent_lengths[d] = n.senescent_length
    # print(n.senescent_lengths)


def mtg_turtle_time_with_constraint(g, day, daily_dynamics, update_visitor=None):  # noqa: ARG001
    """Compute the geometry on each node of the MTG using Turtle geometry.

    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.

    """

    g.properties()["geometry"] = {}
    g.properties()["_plant_translation"] = {}

    max_scale = g.max_scale()

    time = daily_dynamics["Thermal time"]

    distribute_among_organs(g, day, daily_dynamics)
    # print("sen len", g.properties()["senescent_length"].values())

    # to remove
    # visible_leaf_area = g.property("visible_leaf_area")
    # surface1= sum(visible_leaf_area.values())
    

    cereal_visitor = CerealsVisitorGrowth(False)

    # for id in g.vertices():
    #     print(g[id])
    
    
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
            except:  # noqa: E722
                pass
            if g.edge_type(v) == "+":
                turtle.push()
            return True

        def pop_turtle(v):  # noqa: RET503
            n = g.node(v)
            try:
                start_tt = n.start_tt
                if start_tt > time:
                    return False
            except:  # noqa: E722
                pass
            if g.edge_type(v) == "+":
                turtle.pop()

        if g.node(vid).label.startswith("Leaf") or g.node(vid).label.startswith("Stem"):  # noqa: SIM102
            if g.node(vid).start_tt <= time:
                visitor(g, vid, turtle) #, time, growth)
                # turtle.push()
        # plant_id = g.complex_at_scale(vid, scale=1)

        for v in pre_order2_with_filter(g, vid, None, push_turtle, pop_turtle):
            if v == vid:
                continue
            # Done for the leaves
            if g.node(v).start_tt > time:
                continue
            visitor(g, v, turtle) #, time, growth)

        # scene = turtle.getScene()
        return g


    for plant_id in g.component_roots_at_scale_iter(g.root, scale=max_scale):
        g = traverse_with_turtle_time(g, plant_id, time)

    # to remove
    # visible_leaf_area2 = g.property("visible_leaf_area")
    # surface2= sum(visible_leaf_area2.values())
    # dif = daily_dynamics["Plant leaf area"] - surface2
    # real_dif = abs(daily_dynamics["Leaf area increment"] - surface2 + surface1)
    # if real_dif > 1:
    #     print(f'DIFF : {dif} at {time}')
    #     print(surface2-surface1)
    #     print(f'Plant Leaf Area {surface2}')
    #     print(f'Plant Leaf Area STICS {daily_dynamics["Plant leaf area"] }')
    #     print(f'{daily_dynamics["Leaf area increment"]}')

    return g



class CerealsVisitorGrowth(CerealsVisitor):
    def __init__(self, classic):
        super().__init__(classic)

    def __call__(self, g, v, turtle): #, time, growth):
        
        # 1. retrieve the node
        geoms_senesc = g.property("geometry_senescent")
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

        # Try to build the mesh from GC rather than Cylinder
        if n.label.startswith("Leaf") or n.label.startswith("Stem"):
            # update geometry of elements
            # if n.length > 0:
            # print(v)
            # nb_of_growing_internodes, nb_of_growing_leaves = growing_organs(g, time)
            # growth = compute_organ_growth(v, n, time, growth, rate=True)
            mesh = compute_growing_organ_geometry(n, self.classic)
            if mesh:  # To DO : reset to None if calculated so ?
                n.geometry = turtle.transform(mesh)
                n.anchor_point = turtle.getPosition()
                n.heading = turtle.getHeading()
            if v in geoms_senesc and geoms_senesc[v] is not None:
                geoms_senesc[v] = turtle.transform(geoms_senesc[v])

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

        # return growth

