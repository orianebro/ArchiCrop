import matplotlib.pyplot as plt

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
# print(stics_output_data)

file_xml = 'proto_sorghum_plt.xml'
params = ['durvieF', 'ratiodurvieI']

print(read_xml_file(file_xml, params))



height=height_stics[-1]
Smax=LA_stics[-1]
nb_phy=14
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
phyllochron=40
plastochron=phyllochron
leaf_duration=time[-1]/phyllochron-nb_phy # !!!!!!!!!!!!!!!!!!!
stem_duration=leaf_duration
leaf_lifespan=500

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
                    leaf_lifespan)
sorghum.generate_potential_plant()
sorghum.define_development()
growing_plant = sorghum.grow_plant(stics_output_data)

# for val in growing_plant[time[-1]].properties()["senescent_lengths"].values():
#     print([round(v,3) for v in val])

# Pylab plot
# organ dynamics
def plot2d():
    for l in range(len(growing_plant[time[-1]].properties()["leaf_lengths"].values())):
        # print(list(growing_plant[time[-1]].properties()["senescent_lengths"].values())[l])
        # plt.plot(time[1:], [list(growing_plant[t].properties()["leaf_lengths"].values())[l][-1] for t in time[1:]])
        plt.plot(time[1:], [list(growing_plant[t].properties()["leaf_lengths"].values())[l][-1] - list(growing_plant[t].properties()["senescent_lengths"].values())[l][-1] for t in time[1:]])
        plt.plot(time[1:], [la  for la, sen in zip(LA_stics[1:], sen_LA_stics[1:])])
        plt.plot(time[1:], [sen  for la, sen in zip(LA_stics[1:], sen_LA_stics[1:])])
        plt.plot(time[1:], [la - sen for la, sen in zip(LA_stics[1:], sen_LA_stics[1:])])
        
    plt.xlabel("Thermal time")
    plt.ylabel("Leaf length")
    # plt.legend(loc="upper left")

    plt.show()

def plot3d(growing_plant=growing_plant):
    tt = len(time)//2
    scene, _ = build_scene(growing_plant[time[tt]], senescence=True)
    display_scene(scene)

#plot2d()
plot3d()

# print(growing_plant[time[-1]].property_names())
# print(growing_plant[time[-1]].properties()['geometry'])