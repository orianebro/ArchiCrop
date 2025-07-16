from __future__ import annotations

from openalea.archicrop.growth import demand_dist, distribute_to_potential

growing_organs = {
    1:{"potential":10, "visible":5},
    2:{"potential":20, "visible":10},
    3:{"potential":30, "visible":10}
}

incr_0 = distribute_to_potential(growing_organs, 0, demand_dist)
assert(incr_0 == {1: 0.0, 2: 0.0, 3: 0.0})

growing_organs = {
    1:{"potential":10, "visible":5},
    2:{"potential":20, "visible":10},
    3:{"potential":30, "visible":10}
}

incr_50 = distribute_to_potential(growing_organs, 50, demand_dist)
assert(incr_50 == {1: 5.0, 2: 10.0, 3: 20.0})

growing_organs = {
    1:{"potential":10, "visible":5},
    2:{"potential":20, "visible":10},
    3:{"potential":30, "visible":10}
}

incr_35 = distribute_to_potential(growing_organs, 35, demand_dist)
assert(incr_35 == {1: 5.0, 2: 10.0, 3: 20.0})
