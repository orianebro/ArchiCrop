from __future__ import annotations

import matplotlib.pyplot as plt

from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.simulation import read_sti_file, read_xml_file

# Get output data from crop model
stics_output_file = 'mod_smaize.sti'
file_tec_xml = 'Mais_tec.xml'
params_tec = ['densitesem', 'interrang']
tec_stics = read_xml_file(file_tec_xml, params_tec)
sowing_density = tec_stics['densitesem']
inter_row = 0.4
stics_output_data = read_sti_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

file_xml = 'corn_plt.xml'
params_sen = ['durvieF', 'ratiodurvieI']

sen_stics = read_xml_file(file_xml, params_sen)
lifespan = sen_stics['durvieF']
lifespan_early = sen_stics['ratiodurvieI'] * lifespan

for key, value in stics_output_data.items():
    if value["Phenology"] == 'exponential':
        next_key = key + 1
        if next_key in stics_output_data and stics_output_data[next_key]["Phenology"] == 'repro':
            end_veg = value["Thermal time"]
            break

# Set model parameters
height=max(height_stics)
leaf_area=max(LA_stics)
nb_phy=11
wl=0.12
diam_base=2.5 
diam_top=1.5
insertion_angle=30
scurv=0.7
curvature=130
klig=0.6
swmax=0.55
f1=0.64 
f2=0.92
stem_q=1.1
rmax=0.8
skew=0.0005
phyllotactic_angle=180
phyllotactic_deviation=0

leaf_lifespan = [lifespan_early, lifespan]

nb_tillers=0
tiller_delay=1
reduction_factor=1

# Test different values for phyllochron
phyllochrons_to_test = [30]
phyllochrons = [phy for phy in phyllochrons_to_test if end_veg/phy-nb_phy > 1]

for phy in phyllochrons:
    phyllochron = phy
    # print("Phyllochron :", phyllochron)
    plastochron=phyllochron
    # leaf_duration=time[-1]/phy-nb_phy
    # stem_duration=leaf_duration
    # print("Leaf duration :", leaf_duration)

    # Run ArchiCrop model
    maize = ArchiCrop(height=height, 
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
    maize.generate_potential_plant()
    growing_plant = maize.grow_plant()

    for l in range(len(growing_plant[time[-1]].properties()["leaf_lengths"].values())):  # noqa: E741
        # print(list(growing_plant[time[1]].properties()["leaf_lengths"].values()))
        # print(list(growing_plant[time[-1]].properties()["leaf_lengths"].values()))
        plt.plot(time[1:], [list(growing_plant[t].properties()["leaf_lengths"].values())[l][-1] - list(growing_plant[t].properties()["senescent_lengths"].values())[l][-1] for t in time[1:]], color="green", alpha=min(1,20/phyllochron), label=f"{phyllochron}")
    
plt.xlabel("Thermal time")
plt.ylabel("Leaf length")
# plt.legend(loc="upper left")

plt.show()

# print(growing_plant[time[-1]].property_names())
# print(growing_plant[time[-1]].properties()['geometry'])