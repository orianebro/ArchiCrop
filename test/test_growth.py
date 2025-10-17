from __future__ import annotations

import numpy as np
import pytest

from openalea.archicrop import growth
from openalea.archicrop.cereal_plant import cereal

daily_dynamics = {
    "Thermal time": 29.0,
    "Thermal time increment": 29.0,
    "Height increment": 10.0,
    "Leaf area increment": 10.0,
    "Senescent leaf area increment": 10.0
}

def plant(la = 200, h = 100, phyllochron = 30, stem_duration = 2, leaf_duration = 2):
    g = cereal(
        nb_phy=5,
        phyllochron=phyllochron,
        plastochron=phyllochron,
        stem_duration=stem_duration,
        leaf_duration=leaf_duration,
        leaf_lifespan=100,
        end_juv=50,
        nb_tillers=2,
        tiller_delay=1,
        reduction_factor=1,
        height=h,
        leaf_area=la,
        nb_short_phy=2,
        short_phy_height=2,
        wl=0.1,
        diam_base=2.5,
        diam_top=1.5,
        insertion_angle=35,
        scurv=0.7,
        curvature=120,
        klig=0.6,
        swmax=0.55,
        f1=0.64,
        f2=0.92,
        stem_q=1.0,
        rmax=0.7,
        skew=0.005,
        phyllotactic_angle=137.5,
        phyllotactic_deviation=0,
        tiller_angle=30,
        gravitropism_coefficient=0,
        plant_orientation=0,
        spiral=True,
        classic=False
    )

    start = 0.0
    for vid in g.property("leaf_area"): 
        n = g.node(vid)
        n.start_tt = start
        n.end_tt = start + phyllochron * leaf_duration
        start += phyllochron

    return g  # noqa: RET504

def test_init_visible_variables():
    g = plant()
    g = growth.init_visible_variables(g, daily_dynamics)
    g.properties()["visible_length"][1] = 5.0
    g.properties()["stem_diameter"][1] = 2.0
    growth.init_visible_variables(g, daily_dynamics)
    assert all(v == 0.0 for v in g.properties()["visible_length"].values())
    assert all(v == 0.0 for v in g.properties()["stem_diameter"].values())

def test_equal_dist():
    organs = [1, 2, 3]
    result = growth.equal_dist(9, organs)
    assert all(np.isclose(v, 3.0) for v in result.values())

def test_demand_dist():
    organs = {1: {"potential": 2}, 2: {"potential": 4}}
    result = growth.demand_dist(6, organs)
    assert np.isclose(result[1], 2.0)
    assert np.isclose(result[2], 4.0)

def test_get_growing_and_senescing_organs_potential_visible():
    phyllochron = 30
    g = plant(phyllochron=phyllochron)
    g = growth.init_visible_variables(g, daily_dynamics)
    time = phyllochron - 1
    prev_time = 0.0
    growing_internodes, growing_leaves, senescing_leaves = growth.get_growing_and_senescing_organs_potential_visible(g, time, prev_time)
    assert len(growing_leaves) == 1

def test_distribute_to_potential():
    g = plant(phyllochron=30)
    g = growth.init_visible_variables(g, daily_dynamics)
    growing_organs = {
        1:{"potential":10, "visible":5},
        2:{"potential":20, "visible":10},
        3:{"potential":30, "visible":10}
    }
    day=1
    time = daily_dynamics["Thermal time"]
    thermal_time_incr = daily_dynamics["Thermal time increment"]
    increment_to_distribute = daily_dynamics["Leaf area increment"]
    incr_0 = growth.distribute_to_potential(g, growing_organs, day, time, thermal_time_incr, 0.0, growth.demand_dist, growth.update_cereal_leaf_growth_rate)
    assert(incr_0 == {1: 0.0, 2: 0.0, 3: 0.0})

    growing_organs = {
        1:{"potential":10, "visible":5},
        2:{"potential":20, "visible":10},
        3:{"potential":30, "visible":10}
    }
    incr_50 = growth.distribute_to_potential(g, growing_organs, day, time, thermal_time_incr, 50.0, growth.demand_dist, growth.update_cereal_leaf_growth_rate)
    assert(incr_50 == {1: 5.0, 2: 10.0, 3: 20.0})

    growing_organs = {
        1:{"potential":10, "visible":5},
        2:{"potential":20, "visible":10},
        3:{"potential":30, "visible":10}
    }
    incr_35 = growth.distribute_to_potential(g, growing_organs, day, time, thermal_time_incr, 50.0, growth.demand_dist, growth.update_cereal_leaf_growth_rate)
    assert(incr_35 == {1: 5.0, 2: 10.0, 3: 20.0})

def test_distribute_among_organs():
    phyllochron = 30
    g = plant(phyllochron=phyllochron)

    g = growth.init_visible_variables(g, daily_dynamics)
    time = phyllochron - 1
    prev_time = 0
    vid_first_leaf = [vid for vid in g.vertices(scale=g.max_scale()) if g.node(vid).label.startswith("Leaf")][0]  # noqa: RUF015
    pot_first_leaf_area = g.node(vid_first_leaf).leaf_area 
    dynamics = {
        "Thermal time": time,
        "Thermal time increment": time - prev_time,
        "Height increment": 10.0,
        "Leaf area increment": pot_first_leaf_area,
        "Senescent leaf area increment": 0.0
    }
    result = growth.distribute_among_organs(g, 1, dynamics)
    assert result["LA_for_each_leaf"][vid_first_leaf] == pytest.approx(pot_first_leaf_area, rel=0.01)


def test_update_cereal_leaf_growth_area():
    g = plant()
    g = growth.init_visible_variables(g, daily_dynamics)
    leaf_vids = [vid for vid in g.vertices(scale=g.max_scale()) if g.node(vid).label.startswith("Leaf")]
    n = g.node(leaf_vids[0])
    growth.update_cereal_leaf_growth_area(n, 5, 1)
    assert n.visible_leaf_area > 0
    assert n.visible_length > 0

def test_update_cereal_leaf_growth_rate():
    phyllochron = 30
    leaf_duration = 2
    g = plant(phyllochron=phyllochron, leaf_duration=leaf_duration)
    g = growth.init_visible_variables(g, daily_dynamics)
    leaf_vids = [vid for vid in g.vertices(scale=g.max_scale()) if g.node(vid).label.startswith("Leaf")]
    n = g.node(leaf_vids[0])
    dt = phyllochron * leaf_duration - 1
    LA_incr = n.leaf_area - 1
    growth.update_cereal_leaf_growth_rate(n, LA_incr, 1, dt)
    # The increment should not exceed the potential rate * dt
    max_possible_area = n.potential_growth_rate * dt
    assert n.visible_leaf_area <= max_possible_area
    assert n.visible_length > 0
    # If the leaf is not mature, the area should be positive
    assert n.visible_leaf_area > 0

def test_update_stem_growth_height():
    g = plant()
    g = growth.init_visible_variables(g, daily_dynamics)
    stem_vids = [vid for vid in g.vertices(scale=g.max_scale()) if g.node(vid).label.startswith("Stem")]
    n = g.node(stem_vids[0])
    growth.update_stem_growth_height(n, 5, 1)
    assert n.visible_length > 0

def test_update_stem_growth_rate():
    g = plant()
    g = growth.init_visible_variables(g, daily_dynamics)
    stem_vids = [vid for vid in g.vertices(scale=g.max_scale()) if g.node(vid).label.startswith("Stem")]
    n = g.node(stem_vids[0])
    growth.update_stem_growth_rate(n, 4, 1, 2)
    assert n.visible_length > 0

def test_update_leaf_senescence_area():
    g = plant()
    g = growth.init_visible_variables(g, daily_dynamics)
    leaf_vids = [vid for vid in g.vertices(scale=g.max_scale()) if g.node(vid).label.startswith("Leaf")]
    n = g.node(leaf_vids[0])
    growth.update_leaf_senescence_area(n, 5, 1)
    assert n.senescent_area > 0 or n.dead

