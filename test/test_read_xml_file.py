from __future__ import annotations

from openalea.archicrop.simulation import read_xml_file

file_xml = 'sorgho_tec.xml'
params = ['densitesem', 'interrang']

assert(read_xml_file(file_xml, params) == {'densitesem': 10.0, 'interrang': 0.0})

file_xml = 'proto_sorghum_plt.xml'
params = ['durvieF', 'ratiodurvieI']

assert(read_xml_file(file_xml, params) == {'ratiodurvieI': 0.8, 'durvieF': 240.0})