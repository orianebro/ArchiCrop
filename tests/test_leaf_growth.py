import matplotlib.pyplot as plt

from openalea.plantgl.all import surface

from openalea.archicrop.plant_shape import bell_shaped_dist
from openalea.archicrop.geometry import leaf_area_plant, mesh_area
from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.simulation import read_sti_file
from openalea.archicrop.display import build_scene, display_scene

# Get output data from crop model
stics_output_file = 'mod_ssorghum.sti'
sowing_density = 10
inter_row = 0.4
stics_output_data = read_sti_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

# Set model parameters
height=height_stics[-1]
Smax=LA_stics[-1]
nb_phy=15
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

leaf_lifespan = 500

# Test different values for phyllochron
phyllochrons_to_test = [34,36,38,40,42,44,46]
phyllochrons = [phy for phy in phyllochrons_to_test if time[-1]/phy-nb_phy > 0]

for phy in phyllochrons:
    phyllochron = phy
    print("Phyllochron :", phyllochron)
    plastochron=phyllochron
    leaf_duration=time[-1]/phy-nb_phy
    stem_duration=leaf_duration
    print("Leaf duration :", leaf_duration)

    # Run ArchiCrop model
    sorghum = ArchiCrop(height, 
                        nb_phy,
                        Smax,
                        wl, diam_base, diam_top, 
                        insertion_angle, scurv, curvature, 
                        klig, swmax, f1, f2, 
                        stem_q, rmax, skew,
                        phyllotactic_angle,
                        phyllotactic_deviation,
                        phyllochron, 
                        plastochron, 
                        leaf_lifespan,
                        stics_output_data)
    sorghum.generate_potential_plant()
    sorghum.define_development()
    growing_plant = sorghum.grow_plant(stics_output_data)

    # 
    for l in range(len(growing_plant[time[-1]].properties()["leaf_lengths"].values())):
        print(list(growing_plant[time[1]].properties()["leaf_lengths"].values()))
        print(list(growing_plant[time[-1]].properties()["leaf_lengths"].values()))
        plt.plot(time[1:], [list(growing_plant[t].properties()["leaf_lengths"].values())[l][-1] for t in time[1:]], color="green", alpha=min(1,1/leaf_duration), label=f"{phyllochron}")
    
plt.xlabel("Thermal time")
plt.ylabel("Leaf length")
plt.legend(loc="upper left")

plt.show()

# print(growing_plant[time[-1]].property_names())
# print(growing_plant[time[-1]].properties()['geometry'])