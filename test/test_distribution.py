from openalea.archicrop.growth import distribute_to_potential, distribute_among_organs, demand_dist

growing_organs = {
    1:{"potential":10, "visible":5},
    2:{"potential":20, "visible":10},
    3:{"potential":30, "visible":10}
}

increment_to_distribute = 36

increment_for_each_organ = distribute_to_potential(growing_organs, increment_to_distribute, demand_dist)

print(increment_for_each_organ)