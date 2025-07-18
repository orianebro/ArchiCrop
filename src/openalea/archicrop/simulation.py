from __future__ import annotations

from itertools import product

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import qmc

from .archicrop import ArchiCrop
from .cereal_leaf import growing_leaf_area
from .sky_sources import meteo_day


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
        plas = sample[sampled_params.index("plastochron")] if "plastochron" in sampled_params else archi_params["plastochron"]
        nb = sample[sampled_params.index("nb_phy")] if "nb_phy" in sampled_params else archi_params["nb_phy"]
        # print(end_veg/phi-nb)
        if 1 <= end_veg/plas-nb <= 3 and phi < plas:
            sampled_dict = {
                key: int(value) if key in {"nb_phy", "nb_tillers", "nb_short_phy"} else value
                for key, value in zip(sampled_params, sample)
            }
            param_sets.append({**fixed_params, **sampled_dict})
    
    return param_sets


def params_for_curve_fit(param_sets, curves, error_LA, error_height):
    '''
    Select parameters sets for which the model fits the LAI and the height curves of the crop model, with a given error.

    :param param_sets: list of dicts, sets of parameters for ArchiCrop model
    :param curves: 
    :param error_LA: float, error accepted to the LA curve (in cm²)
    :param error_height: float, error accepted to the height curve (in cm)

    :return: dict of lists of parameters, mtg(t), LA(t), height(t) for curve-fitting simulations
    '''
    thermal_time = [value["Thermal time"] for value in curves.values()]
    LA_constraint = [value["Plant leaf area"] for value in curves.values()]
    # sen_LA_constraint = [value["Plant senescent leaf area"] for value in curves.values()]
    height_constraint = [value["Plant height"] for value in curves.values()]
    
    fit_params = []
    non_fit_params = []

    all_sa = []
    all_h = []
    
    fit = True
    for params in param_sets:
        plant = ArchiCrop(daily_dynamics=curves, **params)
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

        # fonction sepraree filter curve
        # Verify that potential leaf area is greater than constraint leaf area
        sa = []
        h = []
        for t in thermal_time:
            sum_area_temp = 0
            for (vid, s,e) in leaf_starts_ends:
                if s <= t < e:
                    sum_area_temp += (t-s)/(e-s) * g.node(vid).leaf_area
                elif t >= e:
                    sum_area_temp += g.node(vid).leaf_area
            sa.append(sum_area_temp)

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
            all_sa.append(sa)
            all_h.append(h)
            fit_params.append(params)

        '''
        growing_plant = plant.grow_plant()
        growing_plant_mtg = list(growing_plant.values())
    
        LA_archicrop = [sum(la - sen 
                            for la, sen in zip(gp.properties()["visible_leaf_area"].values(), gp.properties()["senescent_area"].values())) 
                            for gp in growing_plant_mtg]
        height_archicrop = [sum([vl 
                                 for vid, vl in gp.properties()["visible_length"].items() 
                                 if gp.node(vid).label.startswith("Stem")]) 
                                 for gp in growing_plant_mtg]
    
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
        '''

    return fit_params, non_fit_params, all_sa, all_h


def simulate_fit_params(fit_params, daily_dynamics):
    """ Simulate the growth of plants using the fit parameters."""
    
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


    
def stics_weather_3d(filename, daily_dynamics):
    """Load the weather data from a file and filter it based on the first and last dates of plant growth."""
    df = meteo_day(filename)

    # Get the first and last dates from daily_dynamics
    first_date = list(daily_dynamics.values())[0]["Date"]
    last_date = list(daily_dynamics.values())[-1]["Date"]

    # Use these dates to filter your DataFrame
    return df[(df.daydate >= pd.to_datetime(first_date)) & (df.daydate <= pd.to_datetime(last_date))]