from __future__ import annotations

from random import *  # noqa: F403

import matplotlib.pyplot as plt
from alinea.caribu.data_samples import data_path

from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.display import build_scene
from openalea.archicrop.light_it import compute_light_inter
from openalea.archicrop.stand import agronomic_plot
from openalea.archicrop.stics_io import get_stics_data
from openalea.plantgl.all import Color3, Material

seed(18)  # noqa: F405

# Retrieve STICS management and senescence parameters
tec_stics, stics_output_data, lifespan, lifespan_early = get_stics_data(
    file_tec_xml='Mais_tec.xml',  # Path to the STICS management XML file
    file_plt_xml='corn_plt.xml',  # Path to the STICS plant XML file
    stics_output_file='mod_smaize.sti'  # Path to the STICS output file
)
sowing_density = tec_stics['densitesem']
inter_row = 70

# Retrieve STICS growth and senescence dynamics
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
sen_LA_stics = [value["Plant senescent leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]
par_stics = [value["Absorbed PAR"] for value in stics_output_data.values()]


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
        'plastochron': 40,
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

# Sky
zenith = str(data_path('zenith.light'))

nice_green=Color3((50,100,0))

nplants, positions, domain, domain_area, unit = agronomic_plot(length=1, width=1, sowing_density=sowing_density, inter_row=inter_row, noise=0.1)

inter_plant = 1/sowing_density/inter_row
x_pattern = inter_plant/2
y_pattern = inter_row/2
pattern = (-x_pattern, -y_pattern, x_pattern, y_pattern)
position = (0,0,0)

scenes = {}
for k,v in growing_plant.items():
    scene, nump = build_scene(v, position, leaf_material=Material(nice_green), stem_material=Material(nice_green), senescence=False)
    scenes[k] = scene

par_caribu = []
for scene in scenes.values():
    par_caribu.append(compute_light_inter(scene, zenith, pattern)[0])

plt.plot(time, par_caribu, alpha=0.5, linestyle='--')  # Plot each curve (optional for visualization)
plt.plot(time, par_stics, color="black", label="STICS")

# Labels and legend
plt.xlabel("Thermal time")
plt.ylabel("% of absorbed PAR")
plt.title("Absorbed PAR: 3D canopy vs. STICS")
plt.legend()
plt.show()

# %gui qt
# %run test_use_archicrop.py