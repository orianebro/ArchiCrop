import matplotlib.pyplot as plt

from openalea.plantgl.all import surface

from openalea.archicrop.plant_shape import bell_shaped_dist
from openalea.archicrop.geometry import leaf_area_plant, mesh_area
from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.simulation import read_sti_file, read_xml_file
from openalea.archicrop.display import build_scene, display_scene

stics_output_file = 'mod_ssorghum.sti'
sowing_density = 10
inter_row = 0.4
stics_output_data = read_sti_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

file_xml = 'proto_sorghum_plt.xml'
params_sen = ['durvieF', 'ratiodurvieI']

sen_stics = read_xml_file(file_xml, params_sen)
lifespan = sen_stics['durvieF']
lifespan_early = sen_stics['ratiodurvieI'] * lifespan
leaf_lifespan = [lifespan_early, lifespan]

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

# Test different values for phyllochron
phyllochrons_to_test = [34,36,38,40,42,44,46]
phyllochrons = [phy for phy in phyllochrons_to_test if time[-1]/phy-nb_phy > 0]

fig, ax = plt.subplots(2)

for phy in phyllochrons:
    phyllochron = phy
    print("Phyllochron :", phyllochron)
    plastochron=phyllochron
    leaf_duration=time[-1]/phy-nb_phy
    stem_duration=leaf_duration
    print("Leaf duration :", leaf_duration)

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

    print("LA in :", LA_stics[-1])
    print("LA out theo:", sum(growing_plant[time[-1]].properties()["visible_leaf_area"].values()))
    print("LA out 3D :", sum([surface(geom) for vid, geom in growing_plant[time[-1]].properties()["geometry"].items() if vid in growing_plant[time[-1]].properties()["visible_leaf_area"]])) 
    print("")
    # LA theo =/= real !!!!

    ax[0].plot(range(1,nb_phy+1), growing_plant[time[-1]].properties()["visible_leaf_area"].values(), color="green", alpha=min(1,1/leaf_duration), label=f"{phyllochron}")
    # ax[0].plot(range(1,nb_phy+1), [surface(geom) for vid, geom in growing_plant[time[-1]].properties()["geometry"].items() if vid in growing_plant[time[-1]].properties()["visible_leaf_area"]], color="green", alpha=min(1,1/leaf_duration))
    ax[1].plot(time, [sum(growing_plant[t].properties()["visible_leaf_area"].values()) for t in time], color="green", alpha=min(1,1/leaf_duration), label=f"{phyllochron}")
    # ax[1].plot(time, [sum([surface(geom) for vid, geom in growing_plant[t].properties()["geometry"].items() if vid in growing_plant[time[-1]].properties()["visible_leaf_area"]]) for t in time], color="green", alpha=min(1,1/leaf_duration))

ax[0].plot(range(1,nb_phy+1), bell_shaped_dist(Smax, nb_phy, rmax, skew), color="orange", alpha=0.5, label="Weibull law")
ax[0].set_xlabel("Leaf rank")
ax[0].set_ylabel("Leaf area")
ax[0].legend(loc="upper left")

ax[1].plot(time, LA_stics, "orange", alpha=0.5, label="STICS dynamics")
ax[1].set_xlabel("Thermal time")
ax[1].set_ylabel("Plant leaf area")
ax[1].legend(loc="upper left")

plt.show()

# print(growing_plant[time[-1]].property_names())
# print(growing_plant[time[-1]].properties()['geometry'])