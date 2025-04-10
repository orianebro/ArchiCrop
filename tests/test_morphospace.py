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

scenes = []

phyllochrons_to_test = [34,38,42,46]

for nb_phy in [10,15,20]:
    phyllochrons = [phy for phy in phyllochrons_to_test if time[-1]/phy-nb_phy > 0]
    for phy in phyllochrons:
        phyllochron = phy
        plastochron=phyllochron
        leaf_duration=time[-1]/phy-nb_phy
        stem_duration=leaf_duration

        for insertion_angle in [10,30,50]:

            sorghum = ArchiCrop(height=max(height_stics), 
                                nb_phy=10,
                                Smax=max(LA_stics),
                                wl=0.12, diam_base=2.5, diam_top=1.5, 
                                insertion_angle=insertion_angle, scurv=0.7, curvature=100, 
                                klig=0.6, swmax=0.55, f1=0.64, f2=0.92, 
                                stem_q=1.1, rmax=0.8, skew=0.0005,
                                phyllotactic_angle=180,
                                phyllotactic_deviation=15,
                                phyllochron=phyllochron, 
                                plastochron=plastochron, 
                                leaf_lifespan=leaf_lifespan,
                                increments=stics_output_data)
            sorghum.generate_potential_plant()
            # print("Potential plant generated")
            sorghum.define_development()
            # print("Development defined")
            growing_plant = sorghum.grow_plant(stics_output_data)

            scene = build_scene(growing_plant[time[-1]])
            # display_scene(scene[0])
            scenes.append(scene[0])


# %gui qt
# %run test_use_archicrop.py