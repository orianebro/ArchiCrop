from __future__ import annotations

import os
import sys
import time as t
from multiprocessing import Pool
import datetime

import pandas as pd
import xarray as xr

sys.path.append('../data')
from archi_dict import archi_maize as archi

from openalea.archicrop.simulation import run_simulations
from openalea.mtg.io import write_mtg


def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

def f(n_samples, seed):
    print(f"Running simulation with seed {seed}")
    # Define the inputs for the simulation
    tec_file_xml='../data/Mais_tec.xml'
    plt_file_xml='../data/corn_plt.xml'
    stics_output_file='../data/mod_smaize.sti'
    weather_file = '../data/climaisj.meteo'
    location = {  
    'longitude': 3.87,
    'latitude': 45,
    'altitude': 800,
    'timezone': 'Europe/Paris'}

    # Run the simulation
    daily_dynamics, param_sets, pot_la, pot_h, realized_la, realized_h, nrj_per_plant, mtgs, filters, sowing_density = run_simulations(
        archi_params=archi, 
        tec_file=tec_file_xml, 
        plant_file=plt_file_xml, 
        dynamics_file=stics_output_file, 
        weather_file=weather_file,
        location=location,
        n_samples=n_samples,
        latin_hypercube=True,
        opt_filter_organ_duration=False,
        opt_filter_pot_growth=False,
        opt_filter_realized_growth=False,
        light_inter=False,
        seed=seed)

    # Prepare the data for xarray Dataset
    daily_dyn = {}
    for key in daily_dynamics[1]:
        daily_dyn[key] = [v[key] for v in daily_dynamics.values()]
    dates = pd.to_datetime(daily_dyn['Date'])


    df_archi = pd.DataFrame.from_dict(param_sets, orient='index')
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
    ds.to_netcdf(f"results_{seed}.nc")

if __name__ == '__main__':

    n = 4
    n_samples = [1]*n
    seeds = list(range(1, n+1))

    with Pool(4) as p:
        # info('main line')
        # print(f"Running with {i} CPU")
        start_time = t.time()
        p.starmap_async(f, [(sample, seed) for sample, seed in zip(n_samples, seeds)]).get()
        # p.map(slow_prime_finder, [i*1000000 for i in range(1, 11)])
        end_time = t.time()
        print(f"Time taken with 4 CPU: {end_time - start_time:.2f} seconds")