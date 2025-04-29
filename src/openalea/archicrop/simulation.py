from __future__ import annotations

import xml.etree.ElementTree as ET
from itertools import product

import numpy as np
from scipy.stats import qmc

from .archicrop import ArchiCrop


def dict_ranges_to_all_possible_combinations(d):
    # Convert single floats to lists for consistency
    values = [
        v if isinstance(v, list) else [v]
        for v in d.values()
    ]
    
    # Generate all combinations
    return list(product(*values))


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
        elif isinstance(value, list) and key not in {"leaf_lifespan", "daily_dynamics"}:  # Range to sample
            l_bounds.append(min(value))
            u_bounds.append(max(value))
            
            sampled_params.append(key)

        elif key in {"leaf_lifespan", "daily_dynamics"}:
            fixed_params[key] = value

    
    # Create a Latin Hypercube sampler
    sampler = qmc.LatinHypercube(d=len(l_bounds))  # d = number of parameters
    
    # Generate normalized samples (0 to 1)
    lhs_samples = sampler.random(n=n_samples)
    
    # Scale samples to parameter bounds
    scaled_samples = qmc.scale(lhs_samples, l_bounds=l_bounds, u_bounds=u_bounds)

    # Retrieve phenology
    for key, value in archi_params["daily_dynamics"].items():
        if value["Phenology"] == 'juvenile':
            next_key = key + 1
            # if next_key in archi_params["daily_dynamics"] and archi_params["daily_dynamics"][next_key]["Phenology"] == 'exponential':
            #     end_juv = value["Thermal time"]

        elif value["Phenology"] == 'exponential':
            next_key = key + 1
            if next_key in archi_params["daily_dynamics"] and archi_params["daily_dynamics"][next_key]["Phenology"] == 'repro':
                end_veg = value["Thermal time"]
                break

    # Create parameter sets
    param_sets = []
    for sample in scaled_samples:
        # Combine fixed parameters with sampled parameters
        phi = sample[sampled_params.index("phyllochron")]
        nb = sample[sampled_params.index("nb_phy")]
        if end_veg/phi-nb >= 1:
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
    sen_LA_stics = [value["Senescent leaf area"] for value in curves.values()]
    height_stics = [value["Plant height"] for value in curves.values()]
    
    fitting_sim = {
            'params': [],
            'mtg': [],
            'LA': [],
            'height': []
    }

    
    for params in param_sets:
        sorghum = ArchiCrop(height = max(height_stics), leaf_area = max(LA_stics), **params)
        sorghum.generate_potential_plant()
        growing_plant = sorghum.grow_plant()
        growing_plant_mtg = list(growing_plant.values())
    
        LA_archicrop = [sum(la - sen for la, sen in zip(gp.properties()["visible_leaf_area"].values(), gp.properties()["senescent_area"].values())) for gp in growing_plant_mtg]
        # print(LA_archicrop)
        height_archicrop = [sum([vl for vid, vl in gp.properties()["visible_length"].items() if gp.node(vid).label.startswith("Stem")]) for gp in growing_plant_mtg]
        # --> properties in MTG at plant scale
        # elongation rate (i.e. list of length / day) as organ scale property --> plot
    
        good = True
        for i,(la,h) in enumerate(zip(LA_archicrop, height_archicrop)):
            LA_theo = LA_stics[i] - sen_LA_stics[i]
            if LA_theo*(1-error_LA) <= la <= LA_theo*(1+error_LA) and height_stics[i]*(1-error_height) <= h <= height_stics[i]*(1+error_height): 
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



def read_xml_file(file_xml, params):
    """
    Parses an XML file and retrieves the values of the specified parameters.

    :param file_xml: Path to the XML file.
    :param params: List of parameter names to extract.
    :return: Dictionary with parameter names as keys and extracted values.
    """
    tree = ET.parse(file_xml)
    root = tree.getroot()
    
    result = {}
    
    # Search for all 'param' and 'colonne' elements in the XML
    for elem in root.findall(".//param") + root.findall(".//colonne"):
        param_name = elem.get("nom")  # Get the name attribute
        if param_name in params:
            result[param_name] = float(elem.text.strip()) if elem.text else None
    
    return result


def read_sti_file(file_sti, density):
    """Reads a STICS mod_s*.sti output file and builds a dictionary.
    
    :param file: str, input file of STICS outputs :
        - tempeff(n) : daily efficient thermal time (°C.day)
        - laimax : canopy max LAI (m2/m2)
        - laisen(n) : senescent LAI (m2/m2)
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
    non_zero_height_encountered = False

    with open(file_sti) as file:  # noqa: PTH123
        # Read the header line to get column names
        header = file.readline().strip().split(";")
        # Strip whitespace from column names
        stripped_header = [col.strip() for col in header if col != 'pla']

        # Initialize empty lists for each selected column in the dictionary
        data_dict = {col.strip(): [] for col in stripped_header}

        # Read the rest of the lines (data rows)
        for line in file:
            # Split the line into columns and select only the relevant part
            values = line.strip().split(";")
            values = values[:4] + values[5:]
            # Convert the values to floats
            row = {col.strip(): float(value) for col, value in zip(stripped_header, values)}
            # Check if the height is 0 and break the loop if true, but only after encountering a non-zero height
            if row["hauteur"] != 0.0:
                non_zero_height_encountered = True
                # Append values to the corresponding lists in the dictionary
                for col in stripped_header:
                    data_dict[col.strip()].append(row[col.strip()])
            if non_zero_height_encountered and (row["hauteur"] == 0.0): # or row["laisen(n)"] == 0.0):
                break 

    # start = 21 # 23
    # end = 140
    # density = 10 # density = 20 plants/m2 = 0.002 plants/cm2

    # Thermal time
    thermal_time = list(np.cumsum([float(i) for i in data_dict["tempeff"]]))
    # thermal_time = list(np.cumsum([float(i) for i in data_dict["tmoy(n)"][:end]]))

    # Green LAI
    plant_leaf_area = [10000*float(i)/density for i in data_dict["laimax"]] # from m2/m2 to cm2/plant
    leaf_area_incr = [plant_leaf_area[0]] + [plant_leaf_area[i+1]-plant_leaf_area[i] for i in range(len(plant_leaf_area[1:]))]

    # Senescent LAI
    sen_leaf_area = [10000*float(i)/density for i in data_dict["laisen(n)"]] # from m2/m2 to cm2/plant
    sen_leaf_area_incr = [sen_leaf_area[0]] + [sen_leaf_area[i+1]-sen_leaf_area[i] for i in range(len(sen_leaf_area[1:]))]

    # Phenology
    emergence = data_dict["ilevs"][-1] - data_dict["jul"][0] # from pseudo julian day (from the beginning of the year) to day from begining of the simulation
    end_juv = data_dict["iamfs"][-1] - data_dict["jul"][0]
    max_lai = data_dict["ilaxs"][-1] - data_dict["jul"][0]

    # Height
    height = [float(i)*100 for i in data_dict["hauteur"]] # from m to cm
    height_incr = [height[0]] + [height[i+1]-height[i] for i in range(len(height[1:]))]

    # Absorbed PAR
    par_abs = [float(i)/(0.95*0.48*float(j)) for i, j in zip(data_dict["raint"], data_dict["trg(n)"])] # to % of light intercepted, in MJ/m^2

    return {
        i+1: {"Thermal time": round(thermal_time[i],4),
            "Phenology": 'germination' if i+1 <= emergence else 'juvenile' if emergence < i+1 <= end_juv else 'exponential' if end_juv < i+1 <= max_lai else 'repro',
            "Plant leaf area": round(plant_leaf_area[i],4), 
            "Leaf area increment": round(leaf_area_incr[i],4), 
            "Senescent leaf area": round(sen_leaf_area[i],4),
            "Senescent leaf area increment": round(sen_leaf_area_incr[i],4),
            "Plant height": round(height[i],4), 
            "Height increment": round(height_incr[i],4), 
            "Absorbed PAR": round(par_abs[i],4)}
        for i in range(len(thermal_time))
    }


