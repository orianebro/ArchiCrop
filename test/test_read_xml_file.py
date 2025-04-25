from openalea.archicrop.simulation import read_xml_file

file_xml = 'sorgho_tec.xml'
params = ['densitesem', 'interrang']

print(read_xml_file(file_xml, params))

file_xml = 'proto_sorghum_plt.xml'
params = ['durvieF', 'ratiodurvieI']

print(read_xml_file(file_xml, params))