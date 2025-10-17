from __future__ import annotations

import matplotlib.pyplot as plt

from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.cereal_axis import bell_shaped_dist
from openalea.archicrop.display import build_scene, display_scene
from openalea.archicrop.stics_io import get_stics_data
from openalea.plantgl.all import Color3, Material, surface

# from openalea.archicrop.archi_params import get_archi_params

# Retrieve STICS management and senescence parameters
sowing_density, stics_output_data, lifespan, lifespan_early = get_stics_data(
    file_tec_xml='../data/02NT18SorgV2D1_tec.xml',  # Path to the STICS management XML file
    file_plt_xml='../data/sorgho_imp_M_v10_plt.xml',  # Path to the STICS plant XML file
    stics_output_file='../data/mod_s02NT18SorgV2D1.sti'  # Path to the STICS output file
)

# Retrieve STICS growth and senescence dynamics
dates = [value["Date"] for value in stics_output_data.values()]
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
sen_LA_stics = [value["Plant senescent leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]


# Set ArchiCrop parameters
archi = {
        'height': 3*max(height_stics),
        'leaf_area': 1.2*max(LA_stics),
        'nb_phy': 10,
        'nb_short_phy': 3,
        'wl': 0.08,
        'diam_base': 2.5,
        'diam_top': 1.5,
        'insertion_angle': 35,
        'scurv': 0.7,
        'curvature': 0,
        'klig': 0.6,
        'swmax': 0.55,
        'f1': 0.64,
        'f2': 0.92,
        'stem_q': 1.1,
        'rmax': 0.8,
        'skew': 0.0005,
        'phyllotactic_angle': 180,
        'phyllotactic_deviation': 10,
        'phyllochron': 30,
        # 'plastochron': 40,
        'leaf_lifespan': [lifespan_early, lifespan],
        'nb_tillers': 0,
        'tiller_delay': 1,
        'tiller_angle': 20,
        'reduction_factor': 1,
}

# Instanciate ArchiCrop object
plant = ArchiCrop(daily_dynamics=stics_output_data, **archi)
# Generate a potential plant
plant.generate_potential_plant()
# Simulate growth and senescence of this plant according to the STICS dynamics
growing_plant = plant.grow_plant()


### Plot the 3D scene
times = [t for i,t in enumerate(time) if i%10==0]
mean_time = sum(times) / len(times)
positions = [ (0, 1*(t-mean_time), 0) for t in times]
nice_green = Color3((50, 100, 0))
scene, _ = build_scene(
    [g for i,g in enumerate(list(growing_plant.values())) if i%10==0], 
    position=positions, 
    senescence=True, 
    leaf_material = Material(nice_green), 
    stem_material=Material(nice_green))
display_scene(scene)


### Plot leaf lengths
for l in range(len(growing_plant[dates[-1]].properties()["leaf_lengths"].values())):  # noqa: E741
    # print(list(growing_plant[time[1]].properties()["leaf_lengths"].values()))
    # print(list(growing_plant[time[-1]].properties()["leaf_lengths"].values()))
    plt.plot(time[1:], 
             [list(growing_plant[t].properties()["leaf_lengths"].values())[l][-1] - list(growing_plant[t].properties()["senescent_lengths"].values())[l][-1] for t in dates[1:]], 
             color="green")

plt.xlabel("Thermal time")
plt.ylabel("Leaf length")
# plt.legend(loc="upper left")

plt.show()


### Plot leaf areas constrained vs. realised
fig, ax = plt.subplots(2)

# only works if nb_tillers = 0
ax[0].plot(range(1,plant.nb_phy+1), 
           growing_plant[dates[-1]].properties()["visible_leaf_area"].values(), 
           color="green", label="Realised") 
ax[1].plot(time, 
           [sum(growing_plant[t].properties()["visible_leaf_area"].values()) - sum(growing_plant[t].properties()["senescent_area"].values()) for t in dates], 
           color="green", label="Realised") 

ax[0].plot(range(1,plant.nb_phy+1), 
           bell_shaped_dist(plant.leaf_area, plant.nb_phy, plant.rmax, plant.skew), 
           color="orange", alpha=0.5, label="Potential")
ax[0].set_xlabel("Leaf rank")
ax[0].set_ylabel("Leaf area (cm²)")
ax[0].legend(loc="upper left")

ax[1].plot(time, [la-sen for la, sen in zip(LA_stics, sen_LA_stics)], "orange", alpha=0.5, label="Crop model")
ax[1].set_xlabel("Thermal time (°C.d)")
ax[1].set_ylabel("Plant leaf area (cm²)")
ax[1].legend(loc="upper left")

fig.tight_layout()
# plt.subplots_adjust(bottom=0.2, top=0.8)

plt.show()

# print("LA in :", max(LA_stics))  # noqa: T201
# print("LA out theo:", sum(growing_plant[plant.end_veg].properties()["visible_leaf_area"].values())) # noqa: T201
# print("LA out 3D :", sum( # noqa: T201
#     [surface(geom) 
#      for vid, geom in growing_plant[plant.end_veg].properties()["geometry"].items() 
#      if vid in growing_plant[plant.end_veg].properties()["visible_leaf_area"]]
#      )) 
# LA theo =/= real !!!

# %gui qt
# %run test_use_archicrop.py