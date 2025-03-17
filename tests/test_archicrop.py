from openalea.archicrop.archicrop import ArchiCrop
from openalea.archicrop.simulation import read_sti_file
from openalea.archicrop.display import build_scene, display_scene

stics_output_file = 'mod_ssorghum.sti'
sowing_density = 10
inter_row = 0.4
stics_output_data = read_sti_file(stics_output_file, sowing_density)
time = [value["Thermal time"] for value in stics_output_data.values()]
LA_stics = [value["Plant leaf area"] for value in stics_output_data.values()]
height_stics = [value["Plant height"] for value in stics_output_data.values()]

sorghum = ArchiCrop(height=max(height_stics), 
                    nb_phy=10,
                    Smax=max(LA_stics),
                    wl=0.12, diam_base=2.5, diam_top=1.5, 
                    insertion_angle=30, scurv=0.7, curvature=130, 
                    klig=0.6, swmax=0.55, f1=0.64, f2=0.92, 
                    stem_q=1.1, rmax=0.8, skew=0.0005,
                    phyllotactic_angle=180,
                    phyllotactic_deviation=0,
                    phyllochron=40, 
                    plastochron=40, 
                    leaf_duration=time[-1]/phyllochron-nb_phy, 
                    stem_duration=time[-1]/phyllochron-nb_phy,
                    leaf_lifespan=240)
sorghum.generate_potential_plant()
# print("Potential plant generated")
sorghum.define_development()
# print("Development defined")
growing_plant = sorghum.grow_plant(stics_output_data)


# scene = build_scene(growing_plant[time[-1]])
# display_scene(scene[0])

# %gui qt
# %run test_use_archicrop.py