import matplotlib.pyplot as plt

from openalea.plantgl.all import surface

from openalea.archicrop.plant_shape import bell_shaped_dist
from openalea.archicrop.geometry import leaf_area_plant, mesh_area
from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.simulation import read_sti_file, read_xml_file
from openalea.archicrop.display import build_scene, display_scene

stics_output_file = 'mod_smaize.sti'
file_tec_xml = 'Mais_tec.xml'
params_tec = ['densitesem', 'interrang']
tec_stics = read_xml_file(file_tec_xml, params_tec)
sowing_density = tec_stics['densitesem']
inter_row = 0.4
stics_output_data = read_sti_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
sen_LA_stics = [value["Senescent leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

file_xml = 'corn_plt.xml'
params_sen = ['durvieF', 'ratiodurvieI']

sen_stics = read_xml_file(file_xml, params_sen)
lifespan = sen_stics['durvieF']
lifespan_early = sen_stics['ratiodurvieI'] * lifespan
leaf_lifespan = [lifespan_early, lifespan]

for key, value in stics_output_data.items():
    if value["Phenology"] == 'exponential':
        next_key = key + 1
        if next_key in stics_output_data and stics_output_data[next_key]["Phenology"] == 'repro':
            end_veg = value["Thermal time"]
            break

height=max(height_stics)
Smax=max(LA_stics)
nb_phy=11
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
phyllotactic_deviation=0

# Test different values for phyllochron
phyllochrons_to_test = [30]
phyllochrons = [phy for phy in phyllochrons_to_test if end_veg/phy-nb_phy > 1]

fig, ax = plt.subplots(2)

for phy in phyllochrons:
    phyllochron = phy
    print("Phyllochron :", phyllochron)
    plastochron=phyllochron+11
    # leaf_duration=time[-1]/phy-nb_phy
    # stem_duration=leaf_duration
    # print("Leaf duration :", leaf_duration)

    maize = ArchiCrop(height, 
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
    maize.generate_potential_plant()
    maize.define_development()
    growing_plant = maize.grow_plant()

    print("LA in :", max(LA_stics))
    print("LA out theo:", sum(growing_plant[maize.end_veg].properties()["visible_leaf_area"].values()))
    print("LA out 3D :", sum([surface(geom) for vid, geom in growing_plant[maize.end_veg].properties()["geometry"].items() if vid in growing_plant[maize.end_veg].properties()["visible_leaf_area"]])) 
    print("")
    # LA theo =/= real !!!!

    ax[0].plot(range(1,nb_phy+1), growing_plant[time[-1]].properties()["visible_leaf_area"].values(), color="green", alpha=min(1,20/phy), label="Realised") # label=f"{phyllochron}")
    # ax[0].plot(range(1,nb_phy+1), [surface(geom) for vid, geom in growing_plant[time[-1]].properties()["geometry"].items() if vid in growing_plant[time[-1]].properties()["visible_leaf_area"]], color="green", alpha=min(1,1/leaf_duration))
    ax[1].plot(time, [sum(growing_plant[t].properties()["visible_leaf_area"].values()) - sum(growing_plant[t].properties()["senescent_area"].values()) for t in time], color="green", alpha=min(1,20/phy), label="Realised") # label=f"{phyllochron}")
    # ax[1].plot(time, [sum([surface(geom) for vid, geom in growing_plant[t].properties()["geometry"].items() if vid in growing_plant[time[-1]].properties()["visible_leaf_area"]]) for t in time], color="green", alpha=min(1,1/leaf_duration))

ax[0].plot(range(1,nb_phy+1), bell_shaped_dist(maize.Smax, nb_phy, rmax, skew), color="orange", alpha=0.5, label="Potential")
ax[0].set_xlabel("Leaf rank")
ax[0].set_ylabel("Leaf area (cm²)")
ax[0].legend(loc="upper left")

ax[1].plot(time, [la-sen for la, sen in zip(LA_stics, sen_LA_stics)], "orange", alpha=0.5, label="Crop model")
# ax[1].plot(time, sen_LA_stics, "orange")
ax[1].set_xlabel("Thermal time (°C.d)")
ax[1].set_ylabel("Plant leaf area (cm²)")
ax[1].legend(loc="upper left")

fig.tight_layout()
# plt.subplots_adjust(bottom=0.2, top=0.8)

plt.show()

# print(growing_plant[time[-1]].property_names())
# print(growing_plant[time[-1]].properties()['geometry'])