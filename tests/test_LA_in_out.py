import matplotlib.pyplot as plt

from openalea.archicrop.plant_shape import bell_shaped_dist
from openalea.archicrop.geometry import leaf_area_plant
from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.simulation import retrieve_stics_dynamics_from_file
from openalea.archicrop.display import build_scene, display_scene

stics_output_file = 'mod_ssorghum.sti'
sowing_density = 10
inter_row = 0.4
stics_output_data = retrieve_stics_dynamics_from_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

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

# Test different durations for organ development
leaf_durations = [1,1.5,2,2.5,3,3.5,4]

fig, ax = plt.subplots(2)

for ld in leaf_durations:
    leaf_duration=ld
    stem_duration=leaf_duration
    phyllochron=time[-1]/(nb_phy+leaf_duration)
    print("Phyllochron :", phyllochron)
    plastochron=phyllochron
    leaf_senescence=1000

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
                        leaf_duration, 
                        stem_duration, 
                        leaf_senescence)
    sorghum.generate_potential_plant()
    sorghum.define_development()
    growing_plant = sorghum.grow_plant(stics_output_data)

    print("LA in :", LA_stics[-1])
    print("LA out theo:", sum(growing_plant[time[-1]].properties()["visible_leaf_area"].values()))
    print("LA out 3D :", leaf_area_plant(growing_plant[time[-1]])) 
    # LA theo =/= real !!!!

    ax[0].plot(range(nb_phy), growing_plant[time[-1]].properties()["visible_leaf_area"].values(), color="green", alpha=1/ld)
    ax[1].plot(time, [sum(growing_plant[t].properties()["visible_leaf_area"].values()) for t in time], color="green", alpha=1/ld)

ax[0].plot(range(nb_phy), bell_shaped_dist(Smax, nb_phy, rmax, skew), color="orange", alpha=0.5)
ax[1].plot(time, LA_stics, "orange", alpha=0.5)
plt.show()