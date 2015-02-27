"""
This script creates 1000 different integers of 8 digits each
Each number will be used as a seed for the random resampling
The list of integers will be stored in a file for later use
"""
import random

import folder_assignments as fa

rand_list = []

for i in range(0, 1000):
	rint = random.randint(11111111, 99999999)
	
	if rint not in rand_list:
		rand_list.append(rint)
	else:
		print "oops"


fileloc = "%s/%s/random_seeds.txt"%(fa.mr_folder, fa.numfolder)
if rand_list:
	with open(fileloc, 'w') as fo:
		for rand in rand_list:
			fo.write("%s\n"%rand)
			
else:
	print "there is no rand_list"

print "%s distinct integers have been copied to file %s"%(len(rand_list), fileloc)
