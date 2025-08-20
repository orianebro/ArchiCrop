from __future__ import annotations

import pytest

from openalea.archicrop.cereal_plant import cereal
from openalea.archicrop.light_it import light_interception


@pytest.fixture
def minimal_mtg():
    # Create a minimal MTG with 1 plant, 2 time steps
    g1 = cereal(
        nb_phy=2, phyllochron=30, plastochron=30, stem_duration=2, leaf_duration=2,
        leaf_lifespan=100, end_juv=50, nb_tillers=0, tiller_delay=1, reduction_factor=1,
        height=10, leaf_area=20, nb_short_phy=1, short_phy_height=1, wl=0.1,
        diam_base=2.5, diam_top=1.5, insertion_angle=35, scurv=0.7, curvature=120,
        klig=0.6, swmax=0.55, f1=0.64, f2=0.92, stem_q=1.0, rmax=0.7, skew=0.005,
        phyllotactic_angle=137.5, phyllotactic_deviation=0, tiller_angle=30,
        gravitropism_coefficient=0, plant_orientation=0, spiral=True, classic=False
    )
    g2 = cereal(
        nb_phy=4, phyllochron=30, plastochron=30, stem_duration=2, leaf_duration=2,
        leaf_lifespan=100, end_juv=50, nb_tillers=0, tiller_delay=1, reduction_factor=1,
        height=10, leaf_area=20, nb_short_phy=1, short_phy_height=1, wl=0.1,
        diam_base=2.5, diam_top=1.5, insertion_angle=35, scurv=0.7, curvature=120,
        klig=0.6, swmax=0.55, f1=0.64, f2=0.92, stem_q=1.0, rmax=0.7, skew=0.005,
        phyllotactic_angle=137.5, phyllotactic_deviation=0, tiller_angle=30,
        gravitropism_coefficient=0, plant_orientation=0, spiral=True, classic=False
    )
    # 3 time steps, 2 plants
    return {0:[g1,g1,g1], 1:[g2,g2,g2]}

def test_light_interception_real(minimal_mtg):
    # Minimal daily_dynamics for 2 time steps
    dates = ["1996-04-21", "1996-04-22", "1996-04-23"]
    daily_dynamics = {
        1: {"Date": dates[0], "Incident PAR": 20},
        2: {"Date": dates[1], "Incident PAR": 22},
        3: {"Date": dates[2], "Incident PAR": 24}
    }
    sowing_density = 10
    location = {'longitude': 3.87, 'latitude': 45, 'altitude': 56, 'timezone': 'Europe/Paris'}

    # You need a real weather file compatible with your pipeline
    # For a real test, provide a valid weather file path here:
    weather_file = "climaisj.meteo"  # Make sure this file exists and is valid

    nrj_per_plant = light_interception(
        weather_file=weather_file,
        daily_dynamics=daily_dynamics,
        sowing_density=sowing_density,
        location=location,
        mtgs=minimal_mtg,
        zenith=True,
        save_scenes=False,
        inter_row=70
    )
    # Check output shape and type
    assert len(nrj_per_plant) == len(minimal_mtg)  # 2 plants
    assert len(nrj_per_plant[0]) == len(dates)  # 3 time steps
    assert nrj_per_plant[1] != []