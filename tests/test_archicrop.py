from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.simulation import read_sti_file
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
leaf_lifespan=200

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


# scene = build_scene(growing_plant[time[-1]])
# display_scene(scene[0])

# %gui qt
# %run test_use_archicrop.py