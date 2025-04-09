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
params_sen = ['durvieF', 'ratiodurvieI']

sen_stics = read_xml_file(file_xml, params_sen)
lifespan = sen_stics['durvieF']
lifespan_early = sen_stics['ratiodurvieI'] * lifespan


for key, value in stics_output_data.items():
    if value["Phenology"] == 'exponential':
        next_key = key + 1
        if next_key in stics_output_data and stics_output_data[next_key]["Phenology"] == 'repro':
            time_end_veg = value["Thermal time"]
            break


height=max(height_stics)
Smax=max(LA_stics)
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
rmax=0.8
skew=0.0005
phyllotactic_angle=180
phyllotactic_deviation=0
phyllochron=38
plastochron=phyllochron
leaf_duration=time_end_veg/phyllochron-nb_phy 
stem_duration=leaf_duration
leaf_lifespan=[lifespan_early, lifespan]

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
sorghum.define_development(stics_output_data)
growing_plant = sorghum.grow_plant(stics_output_data)

times = time[-20:-10]
mean_time = sum(times) / len(times)

positions = [ (0, 3*(t-mean_time), 0) for t in times]

scene, _ = build_scene(list(growing_plant.values())[-20:-10], position=positions, senescence=True)
display_scene(scene)

# %gui qt
# %run test_use_archicrop.py