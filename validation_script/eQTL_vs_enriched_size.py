"""
get average ratio of:
eQTL sizes expressed in number of genes
vs
corresponding enrichment list sizes expressed in number of genes
"""

#access genelists
#access entichment
#compare genelist for each trait

from scipy import stats
import numpy as np
import scipy as sp
from numpy import array

import folder_assignments as fa
from data_handlers import read_data

def get_genelist(fn):
	"""
	"""
	trait_genelist = []
	data = read_data(fn)
	for line in data:
		if line.startswith("trait:"):
			trait = line[7:].strip()
		if line.startswith("AT"):
			genelist = line.split()
			trait_genelist.append([trait, genelist])
			
	return trait_genelist



def get_enriched(fn):
	"""
	"""
	datadict = {}
	data = read_data(fn)
	for line in data:
		if line.startswith("trait:"):
			trait = line[7:].strip()
		if line.startswith("AT"):
			genelist = line.split()
			
			datadict[trait] = genelist			
			
	return datadict





def compare_genelists(alist, adict):
	"""
	"""
	ratio_list = []
	for info in alist:
		trait = info[0]
		nr_genelist = len(info[1])
		
		if trait in adict:
			nr_enrichedlist = len(adict[trait])
			
			#ratio = float(nr_genelist) / float(nr_enrichedlist)
			ratio = (float(nr_enrichedlist)) / (float(nr_genelist + nr_enrichedlist))
			#ratio = float(nr_genelist)
			
			ratio_list.append(ratio)
	
	return ratio_list



def process_data(fn):
	"""
	"""
	data = read_data(fn)
	sizes = []
	for line in data:
		if line.startswith("trait:"):
			trait = line[7:].strip()
		if line.startswith("eQTL_size:"):
			size = int(line[11:].strip())
			if size != 0:
				sizes.append(size)
			#if size == 74:
				#print trait

	return sizes


def decriptive_statistics(rlist, dataset, cutoff):
	"""
	"""
	s = array(rlist)
	print "%s \t %s"%(dataset, cutoff)
	print "N : %s"%len(rlist)
	print "Mean : {0:8.6f}".format(s.mean())
	print "Median : {0:8.6f}".format(sp.median(s))
	print "Minimum : {0:8.6f}".format(s.min())
	print "Maximum : {0:8.6f}".format(s.max())
	print "Variance : {0:8.6f}".format(s.var())
	print "Std. deviation : {0:8.6f}".format(s.std())
	print "#####################################"



def main():
	"""
	"""


	#-----------------------------------------------------------
	exp_list = ['Ligterink_2014','Ligterink_2014_gxe','Keurentjes_2007','Snoek_2012']
	#exp_list = ['Ligterink_2014', 'Keurentjes_2007']
	cutoff_list = [3]#,4.3,6.7]
	chromo = [1,2,3,4,5]
	#-----------------------------------------------------------
	
	show_ratio = True
	show_size = False
	
	for dataset in exp_list:
		
		for cutoff in cutoff_list:
			
			if show_ratio:
				############################################################
				####Retrieve genelist
				subfolder_genelist = "genelist_%s"%dataset
				genelist_fn = "%s/%s/%s/genelist_%s_co%s.txt"%(
												fa.mr_folder, fa.gfolder,
												subfolder_genelist, dataset,
												cutoff
												)
				trait_genelist_list = get_genelist(genelist_fn)

				############################################################			
				####Retrieve enriched list
				subfolder_enriched = "enriched_%s"%dataset
				enriched_fn = "%s/%s/%s/enriched_%s_co%s.txt"%(
												fa.mr_folder, fa.enriched_folder,
												subfolder_enriched, dataset,
												cutoff
												)
				dict_trait_enriched = get_enriched(enriched_fn)
				
				############################################################
				####Calculate ratio
				
				ratio_list = compare_genelists(trait_genelist_list, dict_trait_enriched)
				
				
				############################################################
				####Descriptive Statistics
				decriptive_statistics(ratio_list, dataset, cutoff)
				########################################################
				
			if show_size:
				storage_folder = "%s/%s/eQTLsize_%s"%(fa.mr_folder, fa.gfolder, dataset)
				fname = "%s/eQTLsize_%s_co%s.txt"%(storage_folder, dataset, cutoff)
				data_list = process_data(fname)
				decriptive_statistics(data_list, dataset, cutoff)
				############################################################
			
			
			
main()			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
