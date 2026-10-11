from __future__ import annotations

from openalea.archicrop.stics_io import read_xml_file

file_xml = '../data/sorgho_tec.xml'
params = ['densitesem', 'interrang']

assert(read_xml_file(file_xml, params) == {'densitesem': 10.0, 'interrang': 0.0})

file_xml = '../data/plant/proto_sorghum_plt.xml'
params = ['durvieF', 'ratiodurvieI']

assert(read_xml_file(file_xml, params) == {'ratiodurvieI': 0.8, 'durvieF': 240.0})