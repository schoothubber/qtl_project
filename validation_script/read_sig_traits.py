"""
Read traits with a LOD score > 20
dataset is Ligterink_2014
"""

import folder_assignments as fa
from data_handlers import read_data


filename = "%s/%s/lod_of_20.txt"%(fa.mr_folder, fa.raw_folder)
data = read_data(filename)
genelist = []

for line in data:
	if line.startswith("| AT"):
		gene = line[2:12]
		genelist.append(gene)

genelist = sorted(list(set(genelist)))
for g in genelist:
	print g
