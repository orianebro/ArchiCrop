from __future__ import annotations

from random import *  # noqa: F403

import matplotlib.pyplot as plt
from alinea.caribu.data_samples import data_path

from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.display import build_scene
from openalea.archicrop.light_it import compute_light_inter
from openalea.archicrop.simulation import read_sti_file, read_xml_file
from openalea.archicrop.stand import agronomic_plot
from openalea.plantgl.all import Color3, Material

seed(18)  # noqa: F405

file_tec_xml = '../../simulations_STICS/sorgho_tec.xml'
params_tec = ['densitesem', 'interrang']
tec_stics = read_xml_file(file_tec_xml, params_tec)
sowing_density = tec_stics['densitesem']
inter_row = 0.4

stics_output_file = '../../simulations_STICS/mod_ssorghum.sti'
stics_output_data = read_sti_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
sen_LA_stics = [value["Plant senescent leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]
par_stics = [value["Absorbed PAR"] for value in stics_output_data.values()]

file_xml = '../../simulations_STICS/proto_sorghum_plt.xml'
params_sen = ['durvieF', 'ratiodurvieI']

sen_stics = read_xml_file(file_xml, params_sen)
lifespan = sen_stics['durvieF']
lifespan_early = sen_stics['ratiodurvieI'] * lifespan


height=1.5*max(height_stics)
leaf_area=1.5*max(LA_stics)
nb_phy=12
wl=0.12
diam_base=2.5 
diam_top=1.5
insertion_angle=35
scurv=0.7
curvature=100
klig=0.6
swmax=0.55
f1=0.64 
f2=0.92
stem_q=1.1
rmax=0.8
skew=0.0005
phyllotactic_angle=180
phyllotactic_deviation=20
phyllochron=35
plastochron=phyllochron
leaf_lifespan=[lifespan_early, lifespan]
nb_tillers=0
tiller_delay=1
reduction_factor=1

sorghum = ArchiCrop(height=height, 
                    nb_phy=nb_phy,
                    leaf_area=leaf_area,
                    wl=wl, diam_base=diam_base, diam_top=diam_top, 
                    insertion_angle=insertion_angle, scurv=scurv, curvature=curvature, 
                    klig=klig, swmax=swmax, f1=f1, f2=f2, 
                    stem_q=stem_q, rmax=rmax, skew=skew,
                    phyllotactic_angle=phyllotactic_angle,
                    phyllotactic_deviation=phyllotactic_deviation,
                    phyllochron=phyllochron, 
                    plastochron=plastochron, 
                    leaf_lifespan=leaf_lifespan,
                    nb_tillers=nb_tillers, tiller_delay=tiller_delay, reduction_factor=reduction_factor,
                    daily_dynamics=stics_output_data)
sorghum.generate_potential_plant()
growing_plant = sorghum.grow_plant()

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
    par_caribu.append(compute_light_inter(scene, zenith)[0])

plt.plot(time, par_caribu, alpha=0.5, linestyle='--')  # Plot each curve (optional for visualization)
plt.plot(time, par_stics, color="black", label="STICS")
# plt.scatter(time_points, LA_stics)

# Labels and legend
plt.xlabel("Thermal time")
plt.ylabel("% of absorbed PAR")
plt.title("Absorbed PAR: 3D canopy vs. STICS")
plt.legend()
plt.show()

# %gui qt
# %run test_use_archicrop.py