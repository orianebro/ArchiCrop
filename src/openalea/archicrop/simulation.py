from __future__ import annotations

from itertools import product

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


def LHS_param_sampling(archi_params, daily_dynamics, n_samples, seed=42):
    """Generate samples from archi_params dictionary, respecting fixed values."""
    fixed_params = {}
    sampled_params = []

    l_bounds = []
    u_bounds = []

    for key, value in archi_params.items():
        if isinstance(value, (int, float)):  # Fixed parameter
            fixed_params[key] = value
        # Parameter distribution in Latin Hypercube
        elif isinstance(value, list) and key not in {"leaf_lifespan"}:  # Range to sample
            l_bounds.append(min(value))
            u_bounds.append(max(value))
            
            sampled_params.append(key)

        elif key in {"leaf_lifespan"}:
            fixed_params[key] = value

    
    # Create a Latin Hypercube sampler
    sampler = qmc.LatinHypercube(d=len(l_bounds), seed=seed)  # d = number of parameters
    
    # Generate normalized samples (0 to 1)
    lhs_samples = sampler.random(n=n_samples)
    
    # Scale samples to parameter bounds
    scaled_samples = qmc.scale(lhs_samples, l_bounds=l_bounds, u_bounds=u_bounds)

    # Retrieve phenology
    for key, value in daily_dynamics.items():
        if value["Phenology"] == 'juvenile':
            next_key = key + 1
            # if next_key in archi_params["daily_dynamics"] and archi_params["daily_dynamics"][next_key]["Phenology"] == 'exponential':
            #     end_juv = value["Thermal time"]

        elif value["Phenology"] == 'exponential':
            next_key = key + 1
            if next_key in daily_dynamics and daily_dynamics[next_key]["Phenology"] == 'repro':
                end_veg = value["Thermal time"]
                break

    # Create parameter sets
    param_sets = []
    for sample in scaled_samples:
        # Combine fixed parameters with sampled parameters
        phi = sample[sampled_params.index("phyllochron")] if "phyllochron" in sampled_params else archi_params["phyllochron"]
        nb = sample[sampled_params.index("nb_phy")] if "nb_phy" in sampled_params else archi_params["nb_phy"]
        if end_veg/phi-nb >= 1:
            sampled_dict = {
                key: int(value) if key in {"nb_phy", "nb_tillers", "nb_short_phy"} else value
                for key, value in zip(sampled_params, sample)
            }
            param_sets.append({**fixed_params, **sampled_dict})
    
    return param_sets


def regular_param_sampling(archi_params, n_samples):
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
        if end_veg/phi-nb >= 0.5:
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
    sen_LA_stics = [value["Plant senescent leaf area"] for value in curves.values()]
    height_stics = [value["Plant height"] for value in curves.values()]
    
    fitting_sim = {
            'params': [],
            'mtg': [],
            'LA': [],
            'height': []
    }

    non_fitting_sim = {
            'params': [],
            'LA': [],
            'height': []
    }

    
    for params in param_sets:
        plant = ArchiCrop(daily_dynamics=curves, **params)
        plant.generate_potential_plant()
        growing_plant = plant.grow_plant()
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
                non_fitting_sim['params'].append(params)
                non_fitting_sim['LA'].append(LA_archicrop)
                non_fitting_sim['height'].append(height_archicrop)
                break
        if good:
            fitting_sim['params'].append(params)
            fitting_sim['mtg'].append(growing_plant_mtg)
            fitting_sim['LA'].append(LA_archicrop)
            fitting_sim['height'].append(height_archicrop)

    return fitting_sim, non_fitting_sim




