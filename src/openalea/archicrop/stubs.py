import random as rd

from openalea.archicrop.cereals import build_shoot


def generate_potential_plant(archi_parameters):
    potential_plant=build_shoot(archi_parameters)
    return potential_plant

def distribute_among_plants(constraints_crop, nb_of_plants):
    constraint_crop_height=constraints_crop[0]
    constraint_crop_LAI=constraints_crop[1]
    constraint_plants_height=constraint_crop_height
    constraint_plants_LA=constraint_crop_LAI/nb_of_plants
    constraints_plants=[constraint_plants_height, constraint_plants_LA]
    return constraints_plants

def distribute_among_organs(growing_plant, potential_plant, constraints_plants):
    constraint_plants_height=constraints_plants[0]
    constraint_plants_LA=constraints_plants[1] 
    # in growing plant: turtle 
    # identify growing organs(t) = current nb of organs - rank of first growing organ where (n.start_tt <= time < n.end_tt)
    # starting from the lowest (oldest) growing organ:
    # constraint to distribute = constraint_plants_LA
    # constraint for each organ(t) = constraint to distribute / growing organs(t)
    # if organ == leaf: convert constraint for each organ(t) to leaf length
    # if length(t-1) + constraint for each organ(t) <= mature length:
    #   length(t) = length(t-1) + constraint for each organ(t)
    #   constraint to distribute -= constraint for each organ(t)
    # else: 
    #   length(t) = mature length
    #   constraint to distribute -= (mature length - length(t-1))
    #   constraint for each organ(t) = constraint to distribute / growing_organs(t)-1
    growing_plant=42 
    return growing_plant

def grow_plant(potential_plant, constraints_plants):
    plant_at_all_times=[]
    growing_plant=None
    for i in range (len(constraints_plants)):
        growing_plant=distribute_among_organs(growing_plant, potential_plant, constraints_plants)
        plant_at_all_times.append(growing_plant)
    return plant_at_all_times

def spatial_arrangement(crop_parameters):
    crop=42
    return crop

def algo(archi_parameters_ranges, constraints_crop, crop_parameters):
    nb_of_plants=crop_parameters[0]
    plants=[]
    constraints_plants=distribute_among_plants(constraints_crop, nb_of_plants)
    for i in range(nb_of_plants):
        archi_parameters=rd.randrange(archi_parameters_ranges)
        potential_plant=generate_potential_plant(archi_parameters)
        plant_at_all_times=grow_plant(potential_plant, constraints_plants)
        plants.append(plant_at_all_times)
    crop=spatial_arrangement(crop_parameters)
    return crop





