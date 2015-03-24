import sys
import os
import gc

from sortedcontainers import SortedDict

#The following lines enable this script to contact django
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD,Experiment,ExperimentMarker
from qtl.genelist_from_eqtl import (
	marker_logp_list, add_chromosome_and_position, retrieve_logp_above_cutoff, 
	check_region, iterate_marker_tuple, highest_tuple, set_region_from_adjacent_markers, get_chromosome_max, 
	check_regions_on_chromosome, combine_marker_region_data, find_genes_in_regions, normalize_object_data, 
	read_distinct_genes
	)

	


#Main
main_folder = "/home/wouter/xenv_dj1.6/QTL"
script_folder = "validation_script"
result_folder = "validation"
mr_folder =  "%s/%s"%(main_folder, result_folder)
enriched_folder = "enrichment_results"
glist_folder = "gene_lists"
#Raw data
raw_folder = "raw_data"
filename_atreg = "%s/%s/AtRegNet.txt"%(result_folder, raw_folder)
trait_folder = "trait_lists"




def give_genelist_for_trait(trait, cutoff, exp, gxe_boolean):
	"""
	Check enriched genes for each trait for presence of regulator
	
	Perhaps write the genelists to a file for any threshold
	For later use...
	"""
	eQTLsize_list = []
	#print "collecting data for %s of %s with %s"%(trait, exp, cutoff)
	##########################################
	####Run functions towards enrichment######
	##########################################
	
	#Get genelist from eQTL
	l1 = marker_logp_list(trait, gxe_boolean, exp)
	l2 = add_chromosome_and_position(l1)
	l3 = retrieve_logp_above_cutoff(l2, cutoff)
	l4 = check_region(l3)
	l5 = iterate_marker_tuple(l4)
	l6 = set_region_from_adjacent_markers(l5, l2)	
	l7 = combine_marker_region_data(trait, l6)
	gene_dict = find_genes_in_regions(l7)

	#Count the number of found qtls
	qtls = len(gene_dict)
	
	#delete all objects after they become useless
	#because django is a memory horder
	#when gc.collect() is run it will delete all unreferenced dataobjects
	del l1
	del l2
	del l3
	del l4 
	del l5 
	del l6
	del l7
	
	#print "trait %s has %s eQTLs"%(trait, qtls)
	
	if qtls != 0:
		for key, value in gene_dict.iteritems():
			region = "%s_%s - %s_%s"%(key[1], key[2], key[1], key[3])
			eQTLsize_list.append([region, len(value)])
	
	return eQTLsize_list



		

def create_genelist_to_file(fn, dataset, cutoff):
	"""
	"""
	print "creating file for %s"%fn
	with open(fn, 'w') as fo:
		fo.write("dataset: %s"%dataset)
		fo.write("\n")
		fo.write("cutoff: %s"%cutoff)
		fo.write("\n\n")
		

def append_genelist_to_file(fn, gene, eQTL_size_list):
	"""
	
	"""
	print "writing data for %s"%gene
	with open(fn, 'a') as fo:
		fo.write("trait: %s "%gene)
		fo.write("\n")
		for region, eQTL_size in eQTL_size_list:
			fo.write("region: %s\n"%region)
			fo.write("eQTL_size: %s\n"%eQTL_size)
		fo.write("\n")





def get_traits(dataset, cutoff, chromosome):
	"""
	"""
	traitfile = "%s/%s/emr_traitlist_%s_co%s.txt"%(
													mr_folder, 
													trait_folder, 
													dataset, cutoff
													)
	with open(traitfile, 'r') as fo:
		data = fo.readlines()
	
	traitlist = [trait.strip() for trait in data]
	print "%s traits to be processed for %s %s"%(len(traitlist), dataset, cutoff)
	return traitlist

	

def main():
	"""
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/get_eQTL_size.py')		
	"""
	#Datasets
	#exp_list = ['Ligterink_2014','Snoek_2012','Keurentjes_2007']
	exp_list = ['Ligterink_2014']

	#Variables
	chromosome = [1,2,3,4,5]
	cutoff_list = [3, 4.3, 6.7]
	gxe_bool = True

	restart_necessary = False
	
	for dataset in exp_list:
		
		#make a new folder for the dataset
		if not gxe_bool:
			storage_folder = "%s/%s/eQTLsize_%s"%(mr_folder, glist_folder, dataset)
		else:
			storage_folder = "%s/%s/eQTLsize_%s_gxe"%(mr_folder, glist_folder, dataset)
		if not os.path.exists(storage_folder):
			os.mkdir(storage_folder)
		
		for cutoff in cutoff_list:
			
			traitlist = get_traits(dataset, cutoff, chromosome)
			##################################
			if restart_necessary:
				last_trait = "AT3G21420"
				trait_index = traitlist.index(last_trait)
				new_traitlist = traitlist[trait_index:]
			##################################
			print "--------------------------------------------------"
				
			if not gxe_bool:
				fname = "%s/eQTLsize_%s_co%s.txt"%(storage_folder, dataset, cutoff)
			else:
				fname = "%s/eQTLsize_%s_gxe_co%s.txt"%(storage_folder, dataset, cutoff)
			#print fname
			#prep file
			if not os.path.exists(fname):
				create_genelist_to_file(fname, dataset, cutoff)
				print "Created %s"%fname
				print "--------------------------------------------------"
			else:
				print "%s already exists"%fname
				print "--------------------------------------------------"

			
			#get the gene list from the eQTL regions for trait:
			if restart_necessary:
				for trait in new_traitlist:
					#print trait
					eQTL_size_list = give_genelist_for_trait(trait, cutoff, dataset, gxe_bool)
					
					if eQTL_size_list:
						append_genelist_to_file(fname, trait, eQTL_size_list)
						del eQTL_size_list			
						
					gc.collect()
				del traitlist
				del new_traitlist
			else:
				for trait in traitlist:
					#print trait
					eQTL_size_list = give_genelist_for_trait(trait, cutoff, dataset, gxe_bool)
					
					if eQTL_size_list:
						append_genelist_to_file(fname, trait, eQTL_size_list)
						del eQTL_size_list			
						
					gc.collect()
				del traitlist
			
	
	
main()




					
