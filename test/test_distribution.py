from __future__ import annotations

from openalea.archicrop.growth import demand_dist, distribute_to_potential

growing_organs = {
    1:{"potential":10, "visible":5},
    2:{"potential":20, "visible":10},
    3:{"potential":30, "visible":10}
}

increment_to_distribute = 35

increment_for_each_organ = distribute_to_potential(growing_organs, increment_to_distribute, demand_dist)

assert(increment_for_each_organ == {1: 5.0, 2: 10.0, 3: 20.0})