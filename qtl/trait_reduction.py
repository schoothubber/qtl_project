"""
For each dataset get those traits that have a TF in their eQTL
"""

import os
import gc

from qtl.models import Gene,Marker,LOD,Experiment,ExperimentMarker

from qtl.genelist_from_eqtl import (
						marker_logp_list, add_chromosome_and_position, 
						retrieve_logp_above_cutoff, check_region, 
						iterate_marker_tuple, highest_tuple, 
						set_region_from_adjacent_markers, get_chromosome_max, 
						check_regions_on_chromosome, combine_marker_region_data, 
						find_genes_in_regions, normalize_object_data, 
						read_distinct_genes
						)


##################################
#check for each trait if there is a marker higher than cutoff
#if yes than keep the trait




def check_markers(trait, cutoff, exp, gxe_boolean):
	"""
	Check enriched genes for each trait for presence of regulator
	
	Perhaps write the genelists to a file for any threshold
	For later use...
	"""
	
	##########################################
	####Run functions towards enrichment######
	##########################################
	
	#Get genelist from eQTL
	l1 = marker_logp_list(trait, gxe_boolean, exp)
	l2 = add_chromosome_and_position(l1)
	l3 = retrieve_logp_above_cutoff(l2, cutoff)
	
	del l1
	del l2
	
	if l3:
		return l3
		


def get_traitlist(chromosome):
	"""
	"""
	traitlist = []
	for chrom in chromosome:
		traitfile = '/home/wouter/xenv_dj1.6/QTL/validation/trait_lists/traitlist_chr%s.txt'%chrom
		with open(traitfile, 'r') as fo:
			traitdata = fo.readlines()
		
		temp_traitlist = [trait.strip() for trait in traitdata]
		traitlist.extend(temp_traitlist)
	
	return traitlist



def write_trait_to_file(fn, traits):
	"""
	"""
	print "writing %s to file"%fn
	with open(fn, 'w') as fo:
		for trait in traits:
			fo.write(trait)
			fo.write('\n')



def main():
	"""
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/trait_reduction.py')
	"""
	folder = "reduced_traits"
	
	#Datasets
	exp_list = ['Ligterink_2014', 'Ligterink_2014_gxe', 'Keurentjes_2007', 'Snoek_2012']
	cutoff_list = [6.7, 4.3, 3, 2]
	chromosome = [1,2,3,4,5]
	
	###########################################################
	trait_list = get_traitlist(chromosome)
		
	for dataset in exp_list:
		gxe_boolean = False
		
		for cutoff in cutoff_list:
			
			print "Analysing %s with cutoff %s"%(dataset, cutoff)
			
			true_traits = []
			
			subfolder = "%s_%s"%(dataset, cutoff)
			storage = "%s/%s"%(folder, subfolder)
			filename = "%s/%s/true_traits_%s_co%s.txt"%(folder, subfolder, dataset, cutoff)	
			
			if not os.path.exists(storage):
				os.mkdir(storage)	
							
			if dataset == 'Ligterink_2014_gxe':
				dataset = 'Ligterink_2014'
				datset = 'Ligterink_2014_gxe'
				gxe_boolean = True
				filename = "%s/%s/true_traits_%s_co%s.txt"%(folder, subfolder, datset, cutoff)	
				
			for trait in trait_list:
				checks_out_ok = check_markers(trait, cutoff, dataset, gxe_boolean)
				
				if checks_out_ok:
					print trait
					true_traits.append(trait)
					del checks_out_ok
					
			write_trait_to_file(filename, true_traits)
		
		gc.collect()





main()
