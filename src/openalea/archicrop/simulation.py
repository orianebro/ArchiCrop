import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from scipy.stats import qmc

from openalea.plantgl.all import Color3, Material
from oawidgets.plantgl import *
from openalea.mtg.traversal import pre_order2

# from archicrop
from .cereals import build_shoot
from .display import build_scene
from .growth import thermal_time, mtg_turtle_time_with_constraint, init_growth_dict, init_params_for_growth
from .geometry import addSets, leaf_mesh_for_growth, stem_mesh
from .archicrop import ArchiCrop
from .plant_shape import compute_height_growing_plant, compute_leaf_area_growing_plant


def dict_ranges_to_all_possible_combinations(d):
    # Convert single floats to lists for consistency
    values = [
        v if isinstance(v, list) else [v]
        for v in d.values()
    ]
    
    # Generate all combinations
    parameters_combinations = list(product(*values))
    return parameters_combinations


def generate_single_list_dicts(params):
    """
    Generate dictionaries where only one parameter remains a list and all others are single values.
    
    :param params: dict - The input dictionary containing parameters.
    :return: list of dict - A list of dictionaries with only one parameter as a list.
    """
    single_list_dicts = []

    for key, value in params.items():
        if isinstance(value, list):  # Only consider keys with list values
            # Create a base dictionary with single values
            base_dict = {k: (v[0] if isinstance(v, list) else v) for k, v in params.items()}
            # Replace the single value with the original list for the current key
            base_dict[key] = value
            # Add the dictionary to the results
            single_list_dicts.append(base_dict)
    
    return single_list_dicts


def LHS_param_sampling(archi_params, n_samples):
    """Generate samples from archi_params dictionary, respecting fixed values."""
    fixed_params = {}
    sampled_params = []

    l_bounds = []
    u_bounds = []

    for key, value in archi_params.items():
        if isinstance(value, (int, float)):  # Fixed parameter
            fixed_params[key] = value
        # Parameter distribution in Latin Hypercube
        elif isinstance(value, list):  # Range to sample
            l_bounds.append(min(value))
            u_bounds.append(max(value))
            
            sampled_params.append(key)

    
    # Create a Latin Hypercube sampler
    sampler = qmc.LatinHypercube(d=len(l_bounds))  # d = number of parameters
    
    # Generate normalized samples (0 to 1)
    lhs_samples = sampler.random(n=n_samples)
    
    # Scale samples to parameter bounds
    scaled_samples = qmc.scale(lhs_samples, l_bounds=l_bounds, u_bounds=u_bounds)

    # Create parameter sets
    param_sets = []
    for sample in scaled_samples:
        # Combine fixed parameters with sampled parameters
        sampled_dict = {
            key: int(value) if key == "nb_phy" else value
            for key, value in zip(sampled_params, sample)
        }
        param_sets.append({**fixed_params, **sampled_dict})
    
    return param_sets


def params_for_curve_fit(param_sets, curves, error_LA=300, error_height=30):
    '''
    Select parameters sets for which the model fits the LAI and the height curves of the crop model, with a given error.

    :param param_sets: list of dicts, sets of parameters for ArchiCrop model
    :param curves: 
    :param error_LA: float, error accepted to the LA curve (in cm²)
    :param error_height: float, error accepted to the height curve (in cm)

    :return: dict of lists of parameters, mtg(t), LA(t), height(t) for curve-fitting simulations
    '''
    LA_stics = [value["Plant leaf area"] for value in curves.values()]
    height_stics = [value["Plant height"] for value in curves.values()]
    
    fitting_sim = {
            'params': [],
            'mtg': [],
            'LA': [],
            'height': []
    }

    
    for params in param_sets:
        sorghum = ArchiCrop(max(height_stics), **params)
        sorghum.generate_potential_plant()
        sorghum.define_development()
        growing_plant = sorghum.grow_plant(curves)
        growing_plant_mtg = list(growing_plant.values())
    
        LA_archicrop = [sum(compute_leaf_area_growing_plant(gp)) for gp in growing_plant_mtg] 
        height_archicrop = [compute_height_growing_plant(gp) for gp in growing_plant_mtg]
        # --> properties in MTG at plant scale
        # elongation rate (i.e. list of length / day) as organ scale property --> plot
    
        good = True
        for i,(la,h) in enumerate(zip(LA_archicrop, height_archicrop)):
            if LA_stics[i]-error_LA <= la <= LA_stics[i]+error_LA and height_stics[i]-error_height <= h <= height_stics[i]+error_height: 
                good = True
            else:
                good = False
                break
        if good:
            fitting_sim['params'].append(params)
            fitting_sim['mtg'].append(growing_plant_mtg)
            fitting_sim['LA'].append(LA_archicrop)
            fitting_sim['height'].append(height_archicrop)

    return fitting_sim


def generate_potential_plant(nb_phy,
                            height,
                            max_leaf_length,
                            wl,
                            diam_base,
                            diam_top,
                            insertion_angle,
                            scurv,
                            curvature,
                            alpha,
                            stem_q,
                            rmax,
                            skew,
                            phyllotactic_angle,
                            phyllotactic_deviation):
    
    shoot, g = build_shoot(nb_phy,
                            height,
                            max_leaf_length,
                            wl,
                            diam_base,
                            diam_top,
                            insertion_angle,
                            scurv,
                            curvature,
                            alpha,
                            stem_q,
                            rmax,
                            skew,
                            phyllotactic_angle,
                            phyllotactic_deviation)
    return g

def display_plant(g):
    nice_green = Color3((50, 100, 0))
    scene, nump = build_scene(
        g, leaf_material = Material(nice_green), stem_material=Material(nice_green)
    )
    w=PlantGL(scene, group_by_color=True)
    return w 

def generate_and_display_plant(nb_phy,
                                height,
                                max_leaf_length,
                                wl,
                                diam_base,
                                diam_top,
                                insertion_angle,
                                scurv,
                                curvature,
                                alpha,
                                stem_q,
                                rmax,
                                skew,
                                phyllotactic_angle,
                                phyllotactic_deviation):
    shoot, g = build_shoot(nb_phy,
                            height,
                            max_leaf_length,
                            wl,
                            diam_base,
                            diam_top,
                            insertion_angle,
                            scurv,
                            curvature,
                            alpha,
                            stem_q,
                            rmax,
                            skew,
                            phyllotactic_angle,
                            phyllotactic_deviation)
    
    # print("LA plant = ", compute_leaf_area_plant_from_params(nb_phy,
    #                                         max_leaf_length,
    #                                         wl,
    #                                         rmax,
    #                                         skew))

    axes = g.vertices(scale=1)
    metamer_scale = g.max_scale()

    for axis in axes:
        v = next(g.component_roots_at_scale_iter(axis, scale=metamer_scale))

        for metamer in pre_order2(g, v):
            
            n = g.node(metamer)
            n.visible_length = n.mature_length

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
                            inclination=1,
                            stem_diameter=n.stem_diameter,
                        )
                    if n.lrolled > 0:
                        rolled = stem_mesh(
                            n.lrolled, n.lrolled, n.d_rolled, classic=False
                        )
                        if geom is None:
                            geom = rolled
                        else:
                            geom = addSets(rolled, geom, translate=(0, 0, n.lrolled))
            elif n.label.startswith("Stem"):  # stem element
                geom = stem_mesh(n.length, n.visible_length, n.stem_diameter, n.stem_diameter)
    
    nice_green = Color3((50, 100, 0))
    scene, nump = build_scene(
        g, leaf_material = Material(nice_green), stem_material=Material(nice_green), soil_material=Material(Color3((150,100,50)))
    )
    w=PlantGL(scene, group_by_color=True)
    return w 


def retrieve_stics_dynamics_from_file(filename_outputs, density):
    """Reads a STICS mod_s*.sti output file and builds a dictionary.
    
    :param filename_outputs: str, input file of STICS outputs :
        - tmoy(n) : daily thermal time (°C.day)
        - lai(n) : canopy LAI (m2/m2)
        - hauteur : canopy height (m)
        - raint : PAR intercepted (actually, PAR absorbed) by canopy (MJ/m2)
        - trg(n) : global radiation (MJ/m2)
    :return: dict of dicts, for each time step, a dict of values from STICS outputs, converted to be used in ArchiCrop :
        - "Thermal time" (float): thermal time (in °C.day).
        - "Plant leaf area" (float): plant leaf area at a given thermal time (in cm²).
        - "Leaf area increment" (float): leaf area increment at a given thermal time (in cm²).
        - "Plant height" (float): plant height at a given thermal time (in cm).
        - "Height increment" (float): height increment at a given thermal time (in cm).
        - "Absorbed PAR" (float): absorbed PAR at a given thermal time (in MJ/m²)"""
    
    data_dict = {}
    start_column_index = 5
    filter_column = "lai(n)"

    with open(filename_outputs, "r") as file:
        # Read the header line to get column names
        header = file.readline().strip().split(";")
        # Strip whitespace from column names
        stripped_header = [col.strip() for col in header]
        # Select column names starting from the given index
        selected_columns = stripped_header[start_column_index:]
        
        # Check if the filter column exists in the selected columns
        if filter_column not in selected_columns:
            raise ValueError(f"Filter column '{filter_column}' not found in the selected columns.")
        
        # Initialize empty lists for each selected column in the dictionary
        data_dict = {col.strip(): [] for col in selected_columns}
        
        # Read the rest of the lines (data rows)
        for line in file:
            # Split the line into columns and select only the relevant part
            values = line.strip().split(";")[start_column_index:]
            # Convert the values to floats
            row = {col.strip(): float(value) for col, value in zip(selected_columns, values)}
            # Apply the filter condition
            if row[filter_column] != 0:
                # Append values to the corresponding lists in the dictionary
                for col in selected_columns:
                    data_dict[col.strip()].append(row[col.strip()])
    

    # start = 21 # 23
    end = 100 - 23
    # density = 10 # density = 20 plants/m2 = 0.002 plants/cm2

    thermal_time = list(np.cumsum([float(i) for i in data_dict["tempeff"][:end]]))

    leaf_area = [10000*float(i)/density for i in data_dict["lai(n)"][:end]] # from m2/m2 to cm2/plant
    leaf_area_incr = [leaf_area[0]] + [leaf_area[i+1]-leaf_area[i] for i in range(len(leaf_area[1:]))]

    height = [float(i)*100 for i in data_dict["hauteur"][:end]] # from m to cm
    height_incr = [height[0]] + [height[i+1]-height[i] for i in range(len(height[1:]))]

    par_abs = [float(i)/(0.95*0.48*float(j)) for i, j in zip(data_dict["raint"][:end], data_dict["trg(n)"][:end])] # to % of light intercepted, in MJ/m^2

    data = {
        i: {"Thermal time": round(thermal_time[i],4),
            "Plant leaf area": round(leaf_area[i],4), 
            "Leaf area increment": round(leaf_area_incr[i],4), 
            "Plant height": round(height[i],4), 
            "Height increment": round(height_incr[i],4), 
            "Absorbed PAR": round(par_abs[i],4)}
        for i in range(len(thermal_time))
    }

    return data



def model(nb_phy,
            max_leaf_length,
            wl,
            diam_base,
            diam_top,
            insertion_angle,
            scurv,
            curvature,
            alpha,
            stem_q,
            rmax,
            skew,
            phyllotactic_angle,
            phyllotactic_deviation,
            sowing_density, 
            filename_outputs):
    
    nice_green = Color3((50, 100, 0))

    # Retrieve STICS dynamics from files 
    stics_output_data = retrieve_stics_dynamics_from_file(filename_outputs, sowing_density)
    height_potential_plant = max(value["Plant height"] for value in stics_output_data.values())

    # Extract only the increments
    increments = {
        k: {
            "Thermal time": v["Thermal time"],
            "Leaf area increment": v["Leaf area increment"],
            "Height increment": v["Height increment"],
        }
        for k, v in stics_output_data.items()
    }


    # Generate potential plant
    g = generate_potential_plant(nb_phy,
                                height_potential_plant,
                                max_leaf_length,
                                wl,
                                diam_base,
                                diam_top,
                                insertion_angle,
                                scurv,
                                curvature,
                                alpha,
                                stem_q,
                                rmax,
                                skew,
                                phyllotactic_angle,
                                phyllotactic_deviation)

    # Grow population of plants with constraint from crop model
    # and compute light interception for all time steps
    # Add development 
    g = thermal_time(g, phyllochron=50.0, plastochron=40.0, leaf_duration=1.6, stem_duration=1.6)

    # Loop through time
    growing_plant = []
    growth = init_growth_dict()
    g = init_params_for_growth(g)
    for v in increments.values():
        g, growth = mtg_turtle_time_with_constraint(g, v["Thermal time"], v, growth)
        scene, nump = build_scene(
            g, leaf_material = Material(nice_green), stem_material=Material(nice_green)
        )
        growing_plant.append(scene)
    
    # plt.plot(tt_cum[1:], par_caribu, color='orange', label='PAR Caribu')
    # plt.plot(tt_cum[1:], par_stics[1:-1], color='black', label='PAR STICS')
    # plt.xlabel('Thermal time')
    # plt.ylabel('PAR')
    # plt.show()

    return growing_plant


# def grow_pop(g, constraints_crop, sowing_density, inter_row, tt_cum):

#     nice_green = Color3((50, 100, 0))

#     # Add development 
#     g = thermal_time(g)

#     # Spatial arrangement parameters
#     nb_of_plants, positions, domain, domain_area, unit = agronomic_plot(sowing_density=sowing_density, inter_row=inter_row, noise=0.1)

#     scenes = []

#     # Loop through time
#     for i, tt in enumerate(tt_cum):
#         constraints_crop_tt = [constraints_crop[0][i], constraints_crop[1][i]]
#         constraints_plants = distribute_constraint_among_plants(constraints_crop_tt, nb_of_plants)
#         growing_plant = grow_plant_from_constraint(g, tt, constraints_plants)
#         plants = [growing_plant] * nb_of_plants

#         scene, nump = build_scene(plants, positions, leaf_material=Material(nice_green), stem_material=Material(nice_green))
#         # PlantGL(scene, group_by_color=True)
#         scenes.append(scene)
    
#     return scenes



# def grow_pop_and_compute_light_inter(g, constraints_crop, sowing_density, inter_row, tt_cum):

#     nice_green = Color3((50, 100, 0))

#     # Add development 
#     g = thermal_time(g)

#     # Spatial arrangement parameters
#     nb_of_plants, positions, domain, domain_area, unit = agronomic_plot(sowing_density=sowing_density, inter_row=inter_row, noise=0.1)

#     # Loop through time
#     par_caribu = []
#     for i, tt in enumerate(tt_cum):
#         constraints_crop_tt = [constraints_crop[0][i], constraints_crop[1][i]]
#         constraints_plants = distribute_constraint_among_plants(constraints_crop_tt, nb_of_plants)
#         growing_plant = grow_plant_from_constraint(g, tt, constraints_plants)
#         plants = [growing_plant] * nb_of_plants
#         if i > 0:
#             # stand
#             scene, nump = build_scene(plants, positions, leaf_material=Material(nice_green), stem_material=Material(nice_green))
#             # compute light inter
#             par_crop = compute_light_inter(scene)
#             # par_crop = par_crop*0.0145
#             # print(par_crop)
#             par_caribu.append(par_crop)
        
#     return par_caribu


# def model_with_light_inter(nb_phy,
#             max_leaf_length,
#             wl,
#             diam_base,
#             diam_top,
#             insertion_angle,
#             scurv,
#             curvature,
#             alpha,
#             stem_q,
#             rmax,
#             skew,
#             phyllotactic_angle,
#             phyllotactic_deviation,
#             sowing_density, 
#             inter_row,
#             filename_outputs):
    
#     # Retrieve STICS dynamics from files 
#     stics_output_data = retrieve_stics_dynamics_from_file(filename_outputs, sowing_density)
#     height_potential_plant = max(value["Plant height"] for value in stics_output_data.values())


#     # Generate potential plant
#     g = generate_potential_plant(nb_phy,
#                                 height_potential_plant,
#                                 max_leaf_length,
#                                 wl,
#                                 diam_base,
#                                 diam_top,
#                                 insertion_angle,
#                                 scurv,
#                                 curvature,
#                                 alpha,
#                                 stem_q,
#                                 rmax,
#                                 skew,
#                                 phyllotactic_angle,
#                                 phyllotactic_deviation)

#     # Grow population of plants with constraint from crop model
#     # and compute light interception for all time steps
#     par_caribu = grow_pop_and_compute_light_inter(g, 
#                                                   constraints_crop, 
#                                                   sowing_density, 
#                                                   inter_row,
#                                                   tt_cum)
    
#     # plt.plot(tt_cum[1:], par_caribu, color='orange', label='PAR Caribu')
#     # plt.plot(tt_cum[1:], par_stics[1:-1], color='black', label='PAR STICS')
#     # plt.xlabel('Thermal time')
#     # plt.ylabel('PAR')
#     # plt.show()

#     return par_caribu, par_stics