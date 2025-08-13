from __future__ import annotations

import pandas as pd

from openalea.astk.sky_irradiance import sky_irradiance
from openalea.astk.sky_sources import caribu_light_sources, sky_sources

fn = 'climsorj.meteo'
def meteo_day():
    names=['station', 'year', 'month', 'day', 'julian', 'min_temp', 'max_temp', 'rad', 'Penman PET', 'rainfall', 'wind', 'pressure', 'CO2']
    df = pd.read_csv(fn,  header=None, sep='\s+', names=names)  # noqa: PD901
    df["daydate"] = pd.to_datetime(df[["year", "month", "day"]])
    return df

def test_day_to_hour():

    df = meteo_day()  # noqa: PD901
    location ={
    'longitude': 3.87,
    'latitude': 45,
    'altitude': 56,
    'timezone': 'Europe/Paris'}


    for row in df.itertuples():
        irr = sky_irradiance(daydate=row.daydate, day_ghi=row.rad, **location)
        sun, sky = sky_sources(sky_type='blended', sky_irradiance=irr, scale='global')
        lights = caribu_light_sources(sun, sky)
        # then caribu with caribuscene(scene,light=lights,...)
        print(lights)  # noqa: T201


