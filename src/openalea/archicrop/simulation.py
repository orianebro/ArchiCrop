from __future__ import annotations

import os
from datetime import date
from itertools import product

import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import qmc

from openalea.mtg.io import write_mtg

from .archicrop import ArchiCrop
from .cereal_leaf import growing_leaf_area
from .light_it import light_interception
from .stics_io import get_stics_data


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

    archi_params["height"] = 2*max(height_canopy) # [1*max(height_canopy), 2*max(height_canopy)]
    archi_params["leaf_area"] = 2*max(leaf_area_plant) # [1*max(leaf_area_plant), 2*max(leaf_area_plant)]
    archi_params["leaf_lifespan"] = lifespan
    return archi_params


def LHS_param_sampling(archi_params, n_samples, seed=42, latin_hypercube=False):
    """Generate samples from archi_params dictionary, respecting fixed values."""
    fixed_params = {}
    sampled_params = []

    l_bounds = []
    u_bounds = []

    values_lists = []

    for key, value in archi_params.items():
        if isinstance(value, (int, float)):  # Fixed parameter
            fixed_params[key] = value
        # Parameter distribution in Latin Hypercube
        elif isinstance(value, list) and key not in {"leaf_lifespan"}:  # Range to sample

            if latin_hypercube:
                l_bounds.append(min(value))
                u_bounds.append(max(value))
            else:
                # Discretize each parameter
                value_list = np.linspace(min(value), max(value), n_samples)
                values_lists.append(value_list)

            sampled_params.append(key)

        elif key in {"leaf_lifespan"}:
            fixed_params[key] = value


    if latin_hypercube:
        # Create a Latin Hypercube sampler
        sampler = qmc.LatinHypercube(d=len(l_bounds), seed=seed)  # d = number of parameters
        
        # Generate normalized samples (0 to 1)
        lhs_samples = sampler.random(n=n_samples)
        
        # Scale samples to parameter bounds
        samples = qmc.scale(lhs_samples, l_bounds=l_bounds, u_bounds=u_bounds)
    else:
        # Create meshgrid (Cartesian product)
        mesh = np.meshgrid(*values_lists, indexing="ij")

        # Stack into a (N, D) array
        samples = np.stack(mesh, axis=-1).reshape(-1, len(values_lists))

    # Create parameter sets
    param_sets = {}
    for i,sample in enumerate(samples):
        sampled_dict = {
            key: int(value) if key in {"nb_phy", "nb_tillers", "nb_short_phy"} else round(value, 4)
            for key, value in zip(sampled_params, sample)
        }
        param_sets[i] = {**fixed_params, **sampled_dict}
    
    return param_sets


def filter_organ_duration(daily_dynamics, param_sets, opt_filter_organ_duration=True):
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
    # new_param_sets = []
    filters = {}
    for id, sample in param_sets.items():
        if opt_filter_organ_duration:
            phi = sample["phyllochron"] 
            plas = sample["plastochron"] 
            nb = sample["nb_phy"] 
            leaf_duration = end_veg/plas-nb
            stem_duration = end_veg/phi-nb
            # if 0.5 <= leaf_duration <= 3 and 1 <= stem_duration <= 3 and phi < plas:
            if phi < plas:
                # print(leaf_duration)
                filters[id] = {'filter_1' : True}
                # new_param_sets.append(sample)
            else:
                filters[id] = {'filter_1' : False}
        else:
            # filters[id] = {'filter_1' : None}
            filters[id] = {'filter_1' : True}

    return param_sets, filters


def filter_pot_growth(param_sets, daily_dynamics, filters, error_LA, error_height, filter=True):
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

    pot_la = {}
    pot_h = {}
    
    fit = True
    for id, params in param_sets.items():
        if filters[id]['filter_1']:
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
                    fit = False
                    if filter:
                        break

            if filter:
                if fit:
                    pot_la[id] = la
                    pot_h[id] = h
                    filters[id]['filter_2'] = True
                else: 
                    filters[id]['filter_2'] = False
                    pot_la[id] = [None] * len(thermal_time)
                    pot_h[id] = [None] * len(thermal_time)
            else:
                pot_la[id] = la
                pot_h[id] = h
                filters[id]['filter_2'] = None
        else:
            pot_la[id] = [None] * len(thermal_time)
            pot_h[id] = [None] * len(thermal_time)
            filters[id]['filter_2'] = None

    return param_sets, pot_la, pot_h, filters
    

def simulate_fit_params(param_sets, daily_dynamics, filters):
    """ Simulate the growth of plants using the fit parameters.
    Parameters:
        fit_params (list): List of parameter sets.
        daily_dynamics (dict): Dictionary containing daily dynamics data.
    Returns:
        LA_archicrop (list): List of realized leaf area for each parameter set.
        height_archicrop (list): List of realized height for each parameter set.
        mtgs (list): List of MTG objects for each parameter set.
    """
    
    realized_la = {}
    realized_h = {}
    mtgs = {}

    for id, params in param_sets.items():
        if filters[id]['filter_1'] or filters[id]['filter_2']:
            plant = ArchiCrop(daily_dynamics=daily_dynamics, **params)
            plant.generate_potential_plant()
            
            growing_plant = plant.grow_plant()
            growing_plant_mtg = list(growing_plant.values())
            mtgs[id] = growing_plant_mtg

            realized_la[id] = [sum(la - sen 
                                for la, sen in zip(gp.properties()["visible_leaf_area"].values(), gp.properties()["senescent_area"].values())) 
                                for gp in growing_plant_mtg]
            realized_h[id] = [sum([vl 
                                    for vid, vl in gp.properties()["visible_length"].items() 
                                    if gp.node(vid).label.startswith("Stem")]) 
                                    for gp in growing_plant_mtg]
        else:
            realized_la[id] = [None] * len(daily_dynamics)
            realized_h[id] = [None] * len(daily_dynamics)
            mtgs[id] = [None] * len(daily_dynamics)
        
    return realized_la, realized_h, mtgs



def filter_realized_growth(param_sets, realized_la, realized_h, daily_dynamics, mtgs, filters, error_LA=0.05, error_height=0.05, opt_filter_realized_growth=True):
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

    new_realized_la = {}
    new_realized_h = {}
    new_mtgs = {}
    
    for id in param_sets:
        if opt_filter_realized_growth:
            if filters[id]['filter_1'] and filters[id]['filter_2']:
                fit = True
                for i,(la,h) in enumerate(zip(realized_la[id], realized_h[id])):
                    LA_theo = LA_constraint[i] - sen_LA_constraint[i]
                    if LA_theo*(1-error_LA) <= la <= LA_theo*(1+error_LA) and height_constraint[i]*(1-error_height) <= h <= height_constraint[i]*(1+error_height): 
                        fit = True
                    else:
                        fit = False
                        filters[id]['filter_3'] = False
                        new_realized_la[id] = [None] * len(daily_dynamics)
                        new_realized_h[id] = [None] * len(daily_dynamics)
                        new_mtgs[id] = [None] * len(daily_dynamics)
                        break
                if fit:
                    filters[id]['filter_3'] = True
                    new_realized_la[id] = realized_la[id]
                    new_realized_h[id] = realized_h[id]
                    new_mtgs[id] = mtgs[id]
            else:
                filters[id]['filter_3'] = None
                new_realized_la[id] = [None] * len(daily_dynamics)
                new_realized_h[id] = [None] * len(daily_dynamics)
                new_mtgs[id] = [None] * len(daily_dynamics)
        else:
            filters[id]['filter_3'] = None
            new_realized_la[id] = realized_la[id]
            new_realized_h[id] = realized_h[id]
            new_mtgs[id] = mtgs[id]

    return param_sets, new_realized_la, new_realized_h, new_mtgs, filters


def simulate_with_filters(param_sets, daily_dynamics, opt_filter_organ_duration=True, opt_filter_pot_growth=True, opt_filter_realized_growth=True, error_LA_pot=1, error_height_pot=1, error_LA_realized=0.05, error_height_realized=0.05):
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
    param_sets, filters = filter_organ_duration(daily_dynamics=daily_dynamics, param_sets=param_sets, opt_filter_organ_duration=opt_filter_organ_duration)
    param_sets, pot_la, pot_h, filters = filter_pot_growth(param_sets=param_sets, daily_dynamics=daily_dynamics, filters=filters, error_LA=error_LA_pot, error_height=error_height_pot, filter=opt_filter_pot_growth)
    realized_la, realized_h, mtgs = simulate_fit_params(param_sets=param_sets, daily_dynamics=daily_dynamics, filters=filters)
    param_sets, realized_la, realized_h, mtgs, filters = filter_realized_growth(param_sets=param_sets, realized_la=realized_la, realized_h=realized_h, daily_dynamics=daily_dynamics, mtgs=mtgs, filters=filters, error_LA=error_LA_realized, error_height=error_height_realized, opt_filter_realized_growth=opt_filter_realized_growth)
    return param_sets, pot_la, pot_h, realized_la, realized_h, mtgs, filters


def run_simulations(archi_params: dict, 
             tec_file: str, plant_file: str, dynamics_file: str, weather_file: str, location: dict,
             n_samples: int = 100, seed: int = 42, latin_hypercube: bool = False, 
             opt_filter_organ_duration: bool = True, opt_filter_pot_growth: bool = True, opt_filter_realized_growth: bool = True, 
             error_LA_pot: float = 1, error_height_pot: float = 1, error_LA_realized: float = 0.05, error_height_realized: float = 0.05,
             inter_row: float = 70,
             light_inter: bool = True, zenith: bool = False, direct : bool = False, save_scenes: bool = False):

    # Retrieve STICS management and senescence parameters
    sowing_density, daily_dynamics, lifespan, lifespan_early = get_stics_data(
        file_tec_xml=tec_file,  # Path to the STICS management XML file
        file_plt_xml=plant_file,  # Path to the STICS plant XML file
        stics_output_file=dynamics_file  # Path to the STICS output file
    )

    # Complete the architecture parameters with values from daily dynamics.
    archi_params = complete_archi_params(archi_params=archi_params, daily_dynamics=daily_dynamics, lifespan=lifespan, lifespan_early=lifespan_early)
    
    # Sampling parameters using Latin Hypercube Sampling
    param_sets = LHS_param_sampling(archi_params=archi_params, n_samples=n_samples, seed=seed, latin_hypercube=latin_hypercube)

    # Simulate plant growth with fitting parameters
    param_sets, pot_la, pot_h, realized_la, realized_h, mtgs, filters = simulate_with_filters(
        param_sets=param_sets, 
        daily_dynamics=daily_dynamics,
        error_LA_pot=error_LA_pot,
        error_height_pot=error_height_pot, 
        error_LA_realized=error_LA_realized,
        error_height_realized=error_height_realized,
        opt_filter_organ_duration=opt_filter_organ_duration,
        opt_filter_pot_growth=opt_filter_pot_growth,
        opt_filter_realized_growth=opt_filter_realized_growth
    ) 

    # If light_inter is True, compute light interception on plants with fitting parameters
    if light_inter:
        nrj_per_plant = light_interception(
            weather_file=weather_file, 
            daily_dynamics=daily_dynamics, 
            sowing_density=sowing_density, 
            location=location, 
            mtgs=mtgs, 
            zenith=zenith, 
            direct=direct,
            save_scenes=save_scenes, 
            inter_row=inter_row
        )
    else:
        nrj_per_plant = {k : [None] * len(realized_la[k]) for k in realized_la}

    # Dataframe with id, archi_params (dict to df), bool per filter, times series for h, la, nrj + dates ?
    # or xarray

    return daily_dynamics, param_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs, filters, sowing_density

def plot_constained_vs_pot(dates, pot_la, pot_h, leaf_area_plant, height_canopy, sowing_density, stics_color="orange", archicrop_color="green"):
    
    # conversion factor
    cf_cm = 100

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)  # 1 row, 2 columns
    for la in pot_la.values():
        if la[0] is not None:
            axes[0].plot(dates, [a*sowing_density/cf_cm**2 for a in la]) # , color=archicrop_color, alpha=0.6)
    axes[0].plot(dates, [a*sowing_density/cf_cm**2 for a in leaf_area_plant], color=stics_color)
    axes[0].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    axes[0].set_ylabel("LAI (m²/m²)", fontsize=16, fontname="Times New Roman")

    legend_elements_lai = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=2, label='LAI STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='LAI potential morphotypes')
    ]
    axes[0].legend(handles=legend_elements_lai, loc=2, prop={'family': 'Times New Roman', 'size': 12})

    for height in pot_h.values():
        if height[0] is not None:
            axes[1].plot(dates, [h/cf_cm for h in height]) #, color=archicrop_color, alpha=0.6)
    axes[1].plot(dates, [h/cf_cm for h in height_canopy], color=stics_color)
    axes[1].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    axes[1].set_xlabel("Date", fontsize=16, fontname="Times New Roman")
    axes[1].set_ylabel("Crop height (m)", fontsize=16, fontname="Times New Roman")

    legend_elements_height = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=2, label='Height STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='Height potential morphotypes')
    ]
    axes[1].legend(handles=legend_elements_height, loc=2, prop={'family': 'Times New Roman', 'size': 12})

    # Adjust layout
    plt.tight_layout()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_constrained_vs_pot.png")

    # Show the plot
    plt.show()


def plot_constrainted_vs_realized(dates, LA_archicrop, height_archicrop, leaf_area_plant, sen_leaf_area_plant, height_canopy, sowing_density, stics_color="orange", archicrop_color="green"):

    # conversion factor
    cf_cm = 100

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)  # 1 row, 2 columns

    for result in LA_archicrop.values():
        if result[0] is not None:
            axes[0].plot(dates, [r*sowing_density/cf_cm**2 for r in result])  #, color=archicrop_color, alpha=0.6)
    axes[0].plot(dates, [(la-sen)*sowing_density/cf_cm**2 for la, sen in zip(leaf_area_plant, sen_leaf_area_plant)], color=stics_color, alpha=0.9)
    axes[0].set_ylabel("LAI (m²/m²)", fontsize=16, fontname="Times New Roman")
    axes[0].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    # axes[0].set_title("Leaf Area: 3D canopy vs. STICS")
    # axes[0].legend(loc=2)

    legend_elements_lai = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=2, label='LAI STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='LAI morphotypes')
    ]
    axes[0].legend(handles=legend_elements_lai, loc=2, prop={'family': 'Times New Roman', 'size': 12})


    for result in height_archicrop.values():
        if result[0] is not None:
            axes[1].plot(dates, [r/cf_cm for r in result]) #, color=archicrop_color, alpha=0.6)
    axes[1].plot(dates, [h/cf_cm for h in height_canopy], color=stics_color, alpha=0.9)
    axes[1].set_xlabel("Date", fontsize=16, fontname="Times New Roman")
    axes[1].set_ylabel("Crop height (m)", fontsize=16, fontname="Times New Roman")
    axes[0].set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    # axes[1].set_title("Plant height: 3D canopy vs. STICS")

    legend_elements_height = [
        Line2D([0], [0], color=stics_color, alpha=0.9, lw=2, label='Height STICS'),
        Line2D([0], [0], color=archicrop_color, alpha=0.6, lw=2, label='Height morphotypes')
    ]
    axes[1].legend(handles=legend_elements_height, loc=2, prop={'family': 'Times New Roman', 'size': 12})

    # Adjust layout
    plt.tight_layout()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_constrainted_vs_realized.png")

    # Show the plot
    plt.show()


def plot_faPAR(dates, nrj_per_plant, par_incident, par_stics, sowing_density, stics_color="orange", archicrop_color="green"):
    # curves_array = np.array(nrj_per_plant)

    # # Calculate the envelope: min and max values for each time point
    # min_values = curves_array.min(axis=0)
    # max_values = curves_array.max(axis=0)

    # Plotting the envelope along with individual curves for context
    fig, ax = plt.subplots(figsize=(12, 6))
    for curve in nrj_per_plant.values():
        # ax.plot(dates, [nrj*sowing_density/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        ax.plot(dates, [nrj/par for nrj,par in zip(curve, par_incident)]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        # ????????????????

    # ax.fill_between(time_points, min_values, max_values, color="skyblue", alpha=0.4)
    # ax.plot(time_points, min_values, color="blue", linestyle="--", label="Min 3D")
    # ax.plot(time_points, max_values, color="red", linestyle="--", label="Max 3D")
    ax.plot(dates, par_stics, color=stics_color, label="STICS")

    # Labels and legend
    ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Dates") 
    ax.set_ylabel("Fraction of absorbed PAR")
    ax.set_title("Fraction of absorbed PAR: 3D canopy vs. STICS")
    ax.legend()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_faPAR.png")

    plt.show()


def plot_PAR(dates, nrj_per_plant, par_incident, par_stics, sowing_density, stics_color="orange", archicrop_color="green"):
    # curves_array = np.array(nrj_per_plant)

    # # Calculate the envelope: min and max values for each time point
    # min_values = curves_array.min(axis=0)
    # max_values = curves_array.max(axis=0)

    # Plotting the envelope along with individual curves for context
    fig, ax = plt.subplots(figsize=(12, 6))
    for curve in nrj_per_plant.values():
        # ax.plot(dates, [nrj*sowing_density for nrj in curve]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        ax.plot(dates, [nrj for nrj in curve]) #, color=archicrop_color, alpha=0.4, label="ArchiCrop x Caribu")
        # ????????????????

    # ax.fill_between(time_points, min_values, max_values, color="skyblue", alpha=0.4)
    # ax.plot(time_points, min_values, color="blue", linestyle="--", label="Min 3D")
    # ax.plot(time_points, max_values, color="red", linestyle="--", label="Max 3D")
    ax.plot(dates, [abs*inc for abs,inc in zip(par_stics, par_incident)], color=stics_color, label="STICS")

    # Labels and legend
    ax.set_xticks(np.arange(0, len(dates)+1, (len(dates)+1)/8))
    ax.set_xlabel("Dates") 
    ax.set_ylabel("Absorbed PAR")
    ax.set_title("Absorbed PAR: 3D canopy vs. STICS")
    ax.legend()

    # Save figure
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    plt.savefig(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/plot_PAR.png")

    plt.show()


def write_netcdf(daily_dynamics, params_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs, filters, sowing_density, seed):
    # Prepare the data for xarray Dataset
    daily_dyn = {}
    for key in daily_dynamics[1]:
        daily_dyn[key] = [v[key] for v in daily_dynamics.values()]
    dates = pd.to_datetime(daily_dyn['Date'])


    df_archi = pd.DataFrame.from_dict(params_sets, orient='index')
    ds_archi = df_archi.to_xarray().rename({'index':'id'})

    columns_filters = [f"filter_{i}" for i in [1,2,3]]
    df_filters = pd.DataFrame.from_dict(filters, orient='index', columns=columns_filters)
    ds_filters = df_filters.to_xarray().rename({'index':'id'})

    mtgs_string = {k: [write_mtg(g) for g in mtg] for k, mtg in mtgs.items()}

    # Create the xarray Dataset
    ds = xr.Dataset(
        data_vars=dict(  # noqa: C408
            thermal_time = (["time"], daily_dyn["Thermal time"]),
            lai_stics = (["time"], daily_dyn["Plant leaf area"]),
            sen_lai_stics = (["time"], daily_dyn["Plant senescent leaf area"]),
            height_stics = (["time"], daily_dyn["Plant height"]),
            inc_par = (["time"], daily_dyn["Incident PAR"]),
            abs_par_stics = (["time"], daily_dyn["Absorbed PAR"]),
            sowing_density = sowing_density,
            realized_la = (["id", "time"], pd.DataFrame.from_dict(realized_la, orient='index', columns=dates)),
            realized_h = (["id", "time"], pd.DataFrame.from_dict(realized_h, orient='index', columns=dates)),
            pot_la = (["id", "time"], pd.DataFrame.from_dict(pot_la, orient='index', columns=dates)),
            pot_h = (["id", "time"], pd.DataFrame.from_dict(pot_h, orient='index', columns=dates)),
            nrj_per_plant = (["id", "time"], pd.DataFrame.from_dict(nrj_per_plant, orient='index', columns=dates)),
            mtgs = (["id", "time"], pd.DataFrame.from_dict(mtgs_string, orient='index', columns=dates)) 
        ),
        coords=dict(  # noqa: C408
            id = range(len(realized_la)),
            time = dates
        )
    )

    ds = xr.merge([ds, ds_archi, ds_filters])

    # Save the dataset to a NetCDF file
    today_str = date.today().strftime("%Y-%m-%d")
    os.makedirs(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}", exist_ok=True)  # noqa: PTH103
    ds.to_netcdf(f"D:/PhD_Oriane/simulations_ArchiCrop/{today_str}/results_{seed}.nc")