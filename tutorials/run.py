from __future__ import annotations

import sys

import pandas as pd
import xarray as xr

sys.path.append('../data')
from archi_dict import archi_maize as archi

from openalea.archicrop.simulation import run_simulations

# Define the inputs for the simulation
tec_file_xml='Mais_tec.xml'
plt_file_xml='corn_plt.xml'
stics_output_file='mod_smaize.sti'
weather_file = 'climaisj.meteo'
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
    n_samples=10,
    latin_hypercube=True,
    opt_filter_organ_duration=False,
    opt_filter_pot_growth=False,
    opt_filter_realized_growth=False,
    light_inter=False,
    seed=18)

# Prepare the data for xarray Dataset
daily_dyn = {}
for key in daily_dynamics[1]:
    daily_dyn[key] = [v[key] for v in daily_dynamics.values()]
dates = daily_dyn["Date"]

df_archi = pd.DataFrame.from_dict(param_sets, orient='index')
ds_archi = df_archi.to_xarray().rename({'index':'id'})

# Create the xarray Dataset
ds = xr.Dataset(
    data_vars=dict(  # noqa: C408
        # for param in df_archi.columns:
            
        thermal_time = (["time"], daily_dyn["Thermal time"]),
        lai = (["time"], daily_dyn["Plant leaf area"]),
        realized_la = (["id", "time"], pd.DataFrame.from_dict(realized_la, orient='index', columns=dates)),
        realized_h = (["id", "time"], pd.DataFrame.from_dict(realized_h, orient='index', columns=dates)),
    ),
    coords=dict(  # noqa: C408
        id = range(len(realized_la)),
        time = dates
    )
)

ds = xr.merge([ds, ds_archi])


# Save the dataset to a NetCDF file
ds.to_netcdf("../example/simulation_results.nc")

'''
# Read the dataset back from the NetCDF file
ds_read = xr.open_dataset("simulation_results.nc")
# Print the dataset to verify
print(ds_read['lai'].values)
'''
