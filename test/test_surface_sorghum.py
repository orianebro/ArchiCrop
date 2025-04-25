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
sen_LA_stics = [value["Senescent leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

file_xml = 'proto_sorghum_plt.xml'
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
max_LA=max(LA_stics)
nb_phy=12
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
rmax=0.9
skew=0.00005
phyllotactic_angle=180
phyllotactic_deviation=0

# Test different values for phyllochron
phyllochrons_to_test = [20,30,40,50]
phyllochrons = [phy for phy in phyllochrons_to_test if end_veg/phy-nb_phy > 1]

fig, ax = plt.subplots(2)

for phy in phyllochrons:
    phyllochron = phy
    print("Phyllochron :", phyllochron)
    plastochron=phyllochron+10
    # leaf_duration=time[-1]/phy-nb_phy
    # stem_duration=leaf_duration
    # print("Leaf duration :", leaf_duration)

    sorghum = ArchiCrop(height, 
                        nb_phy,
                        max_LA,
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
    # assert(sorghum.Smax == 2*max_LA)
    sorghum.define_development()
    growing_plant = sorghum.grow_plant()

    # times = [t for i,t in enumerate(time) if i%10==0]
    # mean_time = sum(times) / len(times)

    # positions = [ (0, 1*(t-mean_time), 0) for t in times]

    # scene, _ = build_scene([g for i,g in enumerate(list(growing_plant.values())) if i%10==0], position=positions, senescence=True)
    # display_scene(scene)
    # input('next')

#     print("LA in :", LA_stics[-1])
#     print("LA out theo:", sum(growing_plant[sorghum.end_veg].properties()["visible_leaf_area"].values()))
#     print("LA out 3D :", sum([surface(geom) for vid, geom in growing_plant[sorghum.end_veg].properties()["geometry"].items() if vid in growing_plant[sorghum.end_veg].properties()["visible_leaf_area"]])) 
#     print("")
#     # LA theo =/= real !!!!

    ax[0].plot(range(1,nb_phy+1), growing_plant[time[-1]].properties()["visible_leaf_area"].values(), color="green", alpha=min(1,20/phy), label=f"{phyllochron}")
    # ax[0].plot(range(1,nb_phy+1), [surface(geom) for vid, geom in growing_plant[time[-1]].properties()["geometry"].items() if vid in growing_plant[time[-1]].properties()["visible_leaf_area"]], color="green", alpha=min(1,1/leaf_duration))
    ax[1].plot(time, [sum(growing_plant[t].properties()["visible_leaf_area"].values()) - sum(growing_plant[t].properties()["senescent_area"].values()) for t in time], color="green", alpha=min(1,20/phy), label=f"{phyllochron}")
    # ax[1].plot(time, [sum([surface(geom) for vid, geom in growing_plant[t].properties()["geometry"].items() if vid in growing_plant[time[-1]].properties()["visible_leaf_area"]]) for t in time], color="green", alpha=min(1,1/leaf_duration))

ax[0].plot(range(1,nb_phy+1), bell_shaped_dist(sorghum.Smax, nb_phy, rmax, skew), color="orange", alpha=0.5, label="Weibull law")
ax[0].set_xlabel("Leaf rank")
ax[0].set_ylabel("Leaf area")
ax[0].legend(loc="upper left")

ax[1].plot(time, [la-sen for la, sen in zip(LA_stics, sen_LA_stics)], "orange", alpha=0.5, label="STICS dynamics")
ax[1].set_xlabel("Thermal time")
ax[1].set_ylabel("Plant leaf area")
ax[1].legend(loc="upper left")

plt.show()

# print(growing_plant[time[-1]].property_names())
# print(growing_plant[time[-1]].properties()['geometry'])