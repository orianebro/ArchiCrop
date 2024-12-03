import numpy as np
import matplotlib.pyplot as plt
from itertools import product

from openalea.plantgl.all import (Color3, Material, Vector3)
from oawidgets.plantgl import *

# from archicrop
from .cereals import build_shoot
from .display import build_scene
from .stand import agronomic_plot
from .dynamic import thermal_time, grow_plant_from_constraint
from .grow_from_constraint import read_columns_from_file, distribute_constraint_among_plants
from .light_it import compute_light_inter
from .plant_shape import compute_leaf_area_plant_from_params


def dict_ranges_to_all_possible_combinations(d):
    # Convert single floats to lists for consistency
    values = [
        v if isinstance(v, list) else [v]
        for v in d.values()
    ]
    
    # Generate all combinations
    parameters_combinations = list(product(*values))
    return parameters_combinations


def dict_ranges_to_possible_combinations_paired(d):
    
    parameters_combinations = dict_ranges_to_all_possible_combinations(d)

    # Filter combinations where sowing_density == plant_density
    filtered_combinations = [
        combo for combo in parameters_combinations
        if combo[2] == combo[3]  # Ensure sowing_density == plant_density
    ]
    
    return filtered_combinations


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
    
    shoot, potential_plant = build_shoot(nb_phy,
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
    return potential_plant

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
    
    print("LA plant = ", compute_leaf_area_plant_from_params(nb_phy,
                                            max_leaf_length,
                                            wl,
                                            rmax,
                                            skew))
    
    nice_green = Color3((50, 100, 0))
    scene, nump = build_scene(
        g, leaf_material = Material(nice_green), stem_material=Material(nice_green), soil_material=Material(Color3((150,100,50)))
    )
    w=PlantGL(scene, group_by_color=True)
    return w 


def retrieve_stics_dynamics_from_file(filename_outputs, density):
    """STICS outputs : 
        - tmoy(n) : daily thermal time 
        - lai(n) : canopy LAI
        - hauteur : canopy height
        - raint : PAR intercepted (actually, PAR absorbed) by canopy
        - trg(n) : global radiation"""

    columns_outputs = read_columns_from_file(filename_outputs)
    columns_outputs = columns_outputs[5:]

    ## dict with clomn names and then lists (or not), more resilient to changes in order of STICS outputs 

    start = 21 # 21
    end = 90
    # density = 10 # density = 20 plants/m2 = 0.002 plants/cm2

    thermal_time = list(np.cumsum([float(i) for i in columns_outputs[0][start:end-1]]))

    leaf_area = [10000*float(i)/density for i in columns_outputs[1][start:end]] # from m2/m2 to cm2/plant
    leaf_area_incr = [0] + [leaf_area[i+1]-leaf_area[i] for i in range(len(leaf_area[1:]))]

    height = [float(i)*100 for i in columns_outputs[2][start:end]] # from m to cm
    height_incr = [0] + [height[i+1]-height[i] for i in range(len(height[1:]))]

    par_abs = [float(i)/(0.95*0.48*float(j)) for i, j in zip(columns_outputs[3][start:end], columns_outputs[4][start:end])] # to % of light intercepted, in MJ/m^2

    data = {
        thermal_time[i]: {"Plant leaf area": leaf_area[i], 
                          "Leaf area increment": leaf_area_incr[i], 
                          "Plant height": height[i], 
                          "Height increment": height_incr[i], 
                          "Absorbed PAR": par_abs[i]}
        for i in range(len(thermal_time))
    }

    return data


def grow_pop(g, constraints_crop, sowing_density, inter_row, tt_cum):

    nice_green = Color3((50, 100, 0))

    # Add development 
    g = thermal_time(g)

    # Spatial arrangement parameters
    nb_of_plants, positions, domain, domain_area, unit = agronomic_plot(sowing_density=sowing_density, inter_row=inter_row, noise=0.1)

    scenes = []

    # Loop through time
    for i, tt in enumerate(tt_cum):
        constraints_crop_tt = [constraints_crop[0][i], constraints_crop[1][i]]
        constraints_plants = distribute_constraint_among_plants(constraints_crop_tt, nb_of_plants)
        growing_plant = grow_plant_from_constraint(g, tt, constraints_plants)
        plants = [growing_plant] * nb_of_plants

        scene, nump = build_scene(plants, positions, leaf_material=Material(nice_green), stem_material=Material(nice_green))
        # PlantGL(scene, group_by_color=True)
        scenes.append(scene)
    
    return scenes



def grow_pop_and_compute_light_inter(g, constraints_crop, sowing_density, inter_row, tt_cum):

    nice_green = Color3((50, 100, 0))

    # Add development 
    g = thermal_time(g)

    # Spatial arrangement parameters
    nb_of_plants, positions, domain, domain_area, unit = agronomic_plot(sowing_density=sowing_density, inter_row=inter_row, noise=0.1)

    # Loop through time
    par_caribu = []
    for i, tt in enumerate(tt_cum):
        constraints_crop_tt = [constraints_crop[0][i], constraints_crop[1][i]]
        constraints_plants = distribute_constraint_among_plants(constraints_crop_tt, nb_of_plants)
        growing_plant = grow_plant_from_constraint(g, tt, constraints_plants)
        plants = [growing_plant] * nb_of_plants
        if i > 0:
            # stand
            scene, nump = build_scene(plants, positions, leaf_material=Material(nice_green), stem_material=Material(nice_green))
            # compute light inter
            par_crop = compute_light_inter(scene)
            # par_crop = par_crop*0.0145
            # print(par_crop)
            par_caribu.append(par_crop)
        
    return par_caribu


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
            inter_row,
            filename_outputs):
    
    # Retrieve STICS dynamics from files 
    tt_cum, height_cum, la_cum, constraints_crop, par_stics = retrieve_stics_dynamics_from_file(filename_outputs, sowing_density)
    height_potential_plant = max(height_cum)*100


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
    scenes = grow_pop(g, 
                        constraints_crop, 
                        sowing_density,
                        inter_row, 
                        tt_cum)
    
    # plt.plot(tt_cum[1:], par_caribu, color='orange', label='PAR Caribu')
    # plt.plot(tt_cum[1:], par_stics[1:-1], color='black', label='PAR STICS')
    # plt.xlabel('Thermal time')
    # plt.ylabel('PAR')
    # plt.show()

    return scenes



def model_with_light_inter(nb_phy,
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
            inter_row,
            filename_outputs):
    
    # Retrieve STICS dynamics from files 
    tt_cum, height_cum, la_cum, constraints_crop, par_stics = retrieve_stics_dynamics_from_file(filename_outputs, sowing_density)
    height_potential_plant = max(height_cum)*100


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
    par_caribu = grow_pop_and_compute_light_inter(g, 
                                                  constraints_crop, 
                                                  sowing_density, 
                                                  inter_row,
                                                  tt_cum)
    
    # plt.plot(tt_cum[1:], par_caribu, color='orange', label='PAR Caribu')
    # plt.plot(tt_cum[1:], par_stics[1:-1], color='black', label='PAR STICS')
    # plt.xlabel('Thermal time')
    # plt.ylabel('PAR')
    # plt.show()

    return par_caribu, par_stics