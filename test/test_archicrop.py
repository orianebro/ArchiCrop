from __future__ import annotations

from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.display import build_scene, display_scene
from openalea.archicrop.simulation import read_sti_file, read_xml_file
from openalea.plantgl.all import Color3, Material

# Retrieve STICS params about management
# file_tec_xml = 'sorgho_tec.xml'
file_tec_xml = 'Mais_tec.xml'
params_tec = ['densitesem', 'interrang']
tec_stics = read_xml_file(file_tec_xml, params_tec)
sowing_density = tec_stics['densitesem']
inter_row = 0.4

# Retrieve STICS growth and senescence dynamics
# stics_output_file = 'mod_ssorghum.sti'
stics_output_file = 'mod_smaize.sti'
stics_output_data = read_sti_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
sen_LA_stics = [value["Plant senescent leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

# Retrieve STICS params about senescence 
# file_plt_xml = 'proto_sorghum_plt.xml'
file_plt_xml = 'corn_plt.xml'
params_sen = ['durvieF', 'ratiodurvieI']
sen_stics = read_xml_file(file_plt_xml, params_sen)
lifespan = sen_stics['durvieF']
lifespan_early = sen_stics['ratiodurvieI'] * lifespan

# Set ArchiCrop parameters
archi = dict(
    height=3*max(height_stics),
    leaf_area=1.2*max(LA_stics),
    nb_phy=6,
    nb_short_phy=3,
    wl=0.08,
    diam_base=2.5,
    diam_top=1.5,
    insertion_angle=35,
    scurv=0.7,
    curvature=0,
    klig=0.6,
    swmax=0.55,
    f1=0.64,
    f2=0.92,
    stem_q=1.1,
    rmax=0.8,
    skew=0.0005,
    phyllotactic_angle=180,
    phyllotactic_deviation=10,
    phyllochron=30,
    plastochron=40,
    leaf_lifespan=[lifespan_early, lifespan],
    
    nb_tillers=6,
    tiller_delay=1,
    tiller_angle=20,
    reduction_factor=1
)

# Instanciate ArchiCrop object
plant = ArchiCrop(daily_dynamics=stics_output_data, **archi)
# Generate a potential plant
plant.generate_potential_plant()
# Simulate growth and senescence of this plant according to the STICS dynamics
growing_plant = plant.grow_plant()

# Plot the 3D scene
times = [t for i,t in enumerate(time) if i%10==0]
mean_time = sum(times) / len(times)
positions = [ (0, 1*(t-mean_time), 0) for t in times]
nice_green = Color3((50, 100, 0))
scene, _ = build_scene([g for i,g in enumerate(list(growing_plant.values())) if i%10==0], position=positions, senescence=True, leaf_material = Material(nice_green), stem_material=Material(nice_green))
display_scene(scene)

# %gui qt
# %run test_use_archicrop.py