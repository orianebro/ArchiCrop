from __future__ import annotations

from itertools import product

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import qmc

from .archicrop import ArchiCrop
from .cereal_leaf import growing_leaf_area


def leaf_area_plant(g):
    S = 0
    for k,leaf in g.properties()["shape"].items():
        S += growing_leaf_area(leaf, g.properties()["visible_length"][k], g.properties()["mature_length"][k], g.properties()["shape_max_width"][k])
    return S


def dict_ranges_to_all_possible_combinations(d):
    # Convert single floats to lists for consistency
    values = [
        v if isinstance(v, list) else [v]
        for v in d.values()
    ]
    
    # Generate all combinations
    return list(product(*values))


def complete_archi_params(archi_params: dict, daily_dynamics: dict, lifespan: float, lifespan_early: float) -> dict:
    """
    Complete the architecture parameters with values from daily dynamics.
    """
    leaf_area_plant = [value["Plant leaf area"] for value in daily_dynamics.values()]
    height_canopy = [value["Plant height"] for value in daily_dynamics.values()]

    archi_params["height"] = [1.1*max(height_canopy), 5*max(height_canopy)]
    archi_params["leaf_area"] = [1.1*max(leaf_area_plant), 5*max(leaf_area_plant)]
    archi_params["leaf_lifespan"] = [lifespan_early, lifespan]
    return archi_params


def LHS_param_sampling(archi_params, n_samples, seed=42):
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

    # Create parameter sets
    param_sets = []
    for sample in scaled_samples:
        sampled_dict = {
            key: int(value) if key in {"nb_phy", "nb_tillers", "nb_short_phy"} else value
            for key, value in zip(sampled_params, sample)
        }
        param_sets.append({**fixed_params, **sampled_dict})
    
    return param_sets


def filter_organ_duration(daily_dynamics, param_sets):
    """
    Filter the parameter sets to ensure realistic organ duration.
    Parameters:
        daily_dynamics (dict): Dictionary containing daily dynamics data.
        archi_params (dict): Dictionary containing architecture parameters. 
        param_sets (list): List of parameter sets to filter.
    Returns:
        new_param_sets (list): List of parameter sets that meet the organ duration criteria.
    """
    # Get the end of vegetative phase from daily dynamics
    for key, value in daily_dynamics.items():
        if value["Phenology"] == 'exponential':
            next_key = key + 1
            if next_key in daily_dynamics and daily_dynamics[next_key]["Phenology"] == 'repro':
                end_veg = value["Thermal time"]
                break

    # Filter parameter sets to ensure realistic organ duration
    new_param_sets = []
    for sample in param_sets:
        phi = sample["phyllochron"] 
        plas = sample["plastochron"] 
        nb = sample["nb_phy"] 
        leaf_duration = end_veg/plas-nb
        stem_duration = end_veg/phi-nb
        if 1 <= leaf_duration <= 3 and 1 <= stem_duration <= 3 and phi < plas:
            new_param_sets.append(sample)

    return new_param_sets


def filter_pot_growth(param_sets, daily_dynamics, error_LA, error_height):
    '''
    Filter the parameter sets based on potential growth curves.
    Parameters:
        param_sets (list): List of parameter sets to filter.
        daily_dynamics (dict): Dictionary containing daily dynamics data.
        error_LA (float): Acceptable error margin for leaf area.
        error_height (float): Acceptable error margin for height.
    Returns:
        fit_params (list): List of parameter sets that fit the growth curves.
        non_fit_params (list): List of parameter sets that do not fit the growth curves.
        pot_la (list): List of potential leaf area for each parameter set.
        pot_h (list): List of potential height for each parameter set.
    '''
    thermal_time = [value["Thermal time"] for value in daily_dynamics.values()]
    LA_constraint = [value["Plant leaf area"] for value in daily_dynamics.values()]
    # sen_LA_constraint = [value["Plant senescent leaf area"] for value in daily_dynamics.values()]
    height_constraint = [value["Plant height"] for value in daily_dynamics.values()]
    
    fit_params = []
    non_fit_params = []

    pot_la = []
    pot_h = []
    
    fit = True
    for params in param_sets:
        plant = ArchiCrop(daily_dynamics=daily_dynamics, **params)
        plant.generate_potential_plant()
        g = plant.g

        stem_starts_ends = []
        leaf_starts_ends = []
        for vid in g.properties()["start_tt"]:
            if g.node(vid).label.startswith("Stem"):
                stem_starts_ends.append((vid, g.node(vid).start_tt, g.node(vid).end_tt))
            elif g.node(vid).label.startswith("Leaf"):
                leaf_starts_ends.append((vid, g.node(vid).start_tt, g.node(vid).end_tt))
        # Order by start time
        stem_starts_ends.sort(key=lambda x: x[1])  
        leaf_starts_ends.sort(key=lambda x: x[1])

        # Verify that potential leaf area is greater than constraint leaf area
        la = []
        h = []
        for t in thermal_time:
            sum_area_temp = 0
            for (vid, s,e) in leaf_starts_ends:
                if s <= t < e:
                    sum_area_temp += (t-s)/(e-s) * g.node(vid).leaf_area
                elif t >= e:
                    sum_area_temp += g.node(vid).leaf_area
            la.append(sum_area_temp)

            sum_height_temp = 0
            for (vid, s,e) in stem_starts_ends:
                if s <= t < e:
                    sum_height_temp += (t-s)/(e-s) * g.node(vid).mature_length
                elif t >= e:
                    sum_height_temp += g.node(vid).mature_length
            h.append(sum_height_temp)

            if sum_area_temp < LA_constraint[thermal_time.index(t)]*(1-error_LA) or sum_height_temp < height_constraint[thermal_time.index(t)]*(1-error_height):  
                non_fit_params.append(params)
                fit = False
                break

        if fit:
            pot_la.append(la)
            pot_h.append(h)
            fit_params.append(params)

    return fit_params, non_fit_params, pot_la, pot_h
    

def simulate_fit_params(fit_params, daily_dynamics):
    """ Simulate the growth of plants using the fit parameters.
    Parameters:
        fit_params (list): List of parameter sets.
        daily_dynamics (dict): Dictionary containing daily dynamics data.
    Returns:
        LA_archicrop (list): List of realized leaf area for each parameter set.
        height_archicrop (list): List of realized height for each parameter set.
        mtgs (list): List of MTG objects for each parameter set.
    """
    
    LA_archicrop = []
    height_archicrop = []
    mtgs = []

    for params in fit_params:
        plant = ArchiCrop(daily_dynamics=daily_dynamics, **params)
        plant.generate_potential_plant()
        
        growing_plant = plant.grow_plant()
        growing_plant_mtg = list(growing_plant.values())
        mtgs.append(growing_plant_mtg)

        LA_archicrop.append([sum(la - sen 
                            for la, sen in zip(gp.properties()["visible_leaf_area"].values(), gp.properties()["senescent_area"].values())) 
                            for gp in growing_plant_mtg])
        height_archicrop.append([sum([vl 
                                for vid, vl in gp.properties()["visible_length"].items() 
                                if gp.node(vid).label.startswith("Stem")]) 
                                for gp in growing_plant_mtg])
        
    return LA_archicrop, height_archicrop, mtgs



def filter_realized_growth(param_sets, LA_archicrop, height_archicrop, daily_dynamics, error_LA=0.05, error_height=0.05):
    '''
    Filter the parmeters sets based on the deviation of the realized growth from the constraints, within a given error margin.
    Parameters: 
        param_sets (list): List of parameter sets to filter.
        daily_dynamics (dict): Dictionary containing daily dynamics data.
        error_LA (float): Acceptable error margin for leaf area.
        error_height (float): Acceptable error margin for height.
    Returns:
        fit_params (list): List of parameter sets that fit the growth curves.
        non_fit_params (list): List of parameter sets that do not fit the growth curves.    
        realized_la (list): List of realized leaf area for each parameter set.
        realized_h (list): List of realized height for each parameter set.
    '''
    # thermal_time = [value["Thermal time"] for value in daily_dynamics.values()]
    LA_constraint = [value["Plant leaf area"] for value in daily_dynamics.values()]
    sen_LA_constraint = [value["Plant senescent leaf area"] for value in daily_dynamics.values()]
    height_constraint = [value["Plant height"] for value in daily_dynamics.values()]

    fit_params = []
    non_fit_params = []

    realized_la = []
    realized_h = []
    
    for (params, leaf_areas, heights) in zip(param_sets, LA_archicrop, height_archicrop):
    
        fit = True
        for i,(la,h) in enumerate(zip(leaf_areas, heights)):
            LA_theo = LA_constraint[i] - sen_LA_constraint[i]
            if LA_theo*(1-error_LA) <= la <= LA_theo*(1+error_LA) and height_constraint[i]*(1-error_height) <= h <= height_constraint[i]*(1+error_height): 
                fit = True
            else:
                fit = False
                non_fit_params.append(params)
                break
        if fit:
            fit_params.append(params)
            realized_la.append(leaf_areas)
            realized_h.append(heights)

    return fit_params, non_fit_params, realized_la, realized_h


def simulate_with_filters(param_sets, daily_dynamics, opt_filter_pot_growth=True, opt_filter_realized_growth=True, error_LA=0.05, error_height=0.05):
    """
    Simulate the growth of plants using the parameter sets and filter them based on realized growth.
    Parameters:
        param_sets (list): List of parameter sets to simulate.
        daily_dynamics (dict): Dictionary containing daily dynamics data.
        error_LA (float): Acceptable error margin for leaf area.
        error_height (float): Acceptable error margin for height.
    Returns:
        fit_params (list): List of parameter sets that fit the growth curves.
        non_fit_params (list): List of parameter sets that do not fit the growth curves.
        realized_la (list): List of realized leaf area for each parameter set.
        realized_h (list): List of realized height for each parameter set.
    """
    non_fit_params = []
    fit_params = filter_organ_duration(daily_dynamics, param_sets)
    if opt_filter_pot_growth:
        fit_params, non_fit_params_temp, pot_la, pot_h = filter_pot_growth(fit_params, daily_dynamics, error_LA, error_height)
        non_fit_params.extend(non_fit_params_temp)
    realized_la, realized_h, mtgs = simulate_fit_params(fit_params, daily_dynamics)
    if opt_filter_realized_growth:
        fit_params, non_fit_params_temp, realized_la, realized_h = filter_realized_growth(fit_params, realized_la, realized_h, daily_dynamics, error_LA, error_height)
        non_fit_params.extend(non_fit_params_temp)

    return fit_params, non_fit_params, realized_la, realized_h, mtgs


def plot_constrainted_vs_realized(thermal_time, LA_archicrop, height_archicrop, leaf_area_plant, sen_leaf_area_plant, height_canopy, sowing_density):

    # conversion factor
    cf_cm = 100

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)  # 1 row, 2 columns

    for result in LA_archicrop:
        axes[0].plot(thermal_time, [r*sowing_density/cf_cm**2 for r in result],  color="green", alpha=0.6)
    axes[0].plot(thermal_time, [(la-sen)*sowing_density/cf_cm**2 for la, sen in zip(leaf_area_plant, sen_leaf_area_plant)], color="black", alpha=0.9)
    axes[0].set_ylabel("LAI (m²/m²)", fontsize=16, fontname="Times New Roman")
    # axes[0].set_title("Leaf Area: 3D canopy vs. STICS")
    # axes[0].legend(loc=2)

    legend_elements_lai = [
        Line2D([0], [0], color='black', alpha=0.9, lw=2, label='LAI STICS'),
        Line2D([0], [0], color='green', alpha=0.6, lw=2, label='LAI morphotypes')
    ]
    axes[0].legend(handles=legend_elements_lai, loc=2, prop={'family': 'Times New Roman', 'size': 12})


    for result in height_archicrop:
        axes[1].plot(thermal_time, [r*0.01 for r in result], color="orange", alpha=0.6)
    axes[1].plot(thermal_time, [h*0.01 for h in height_canopy], color="black", alpha=0.9)
    axes[1].set_xlabel("Thermal time (°C.day)", fontsize=16, fontname="Times New Roman")
    axes[1].set_ylabel("Crop height (m)", fontsize=16, fontname="Times New Roman")
    # axes[1].set_title("Plant height: 3D canopy vs. STICS")

    legend_elements_height = [
        Line2D([0], [0], color='black', alpha=0.9, lw=2, label='Height STICS'),
        Line2D([0], [0], color='orange', alpha=0.6, lw=2, label='Height morphotypes')
    ]
    axes[1].legend(handles=legend_elements_height, loc=2, prop={'family': 'Times New Roman', 'size': 12})

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()

