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
from qtl.go_enrichment import (
					read_once_csv, segregate_gene_list, annotate_from_csv, 
					read_go_info_once_csv, unique_GO_list, lets_see_tables, 
					flatten_array, count_all_goterms, populate_contingency_table, 
					fish_for_python, extract_significant_result_fish,
					extract_significant_result_mult,
					correct_pvalues_for_multiple_testing, make_qtl_go_dict, 
					post_process_A, post_process_B
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




def give_enrichment_for_trait(trait, cutoff, exp, gxe_boolean, fisher_alpha_name, mult_alpha_name, postA_name, postB_name, data_dict):
	"""
	Check enriched genes for each trait for presence of regulator
	
	Perhaps write the genelists to a file for any threshold
	For later use...
	"""
	
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
	
		gene_list = read_distinct_genes(gene_dict)
		
		absent_genes, present_genes = segregate_gene_list(data_dict, gene_list)	
		gene_inside_qtl_dict, gene_outside_qtl_dict, tot_in_qtl, tot_out_qtl = annotate_from_csv(data_dict, present_genes)

		golist_unique = unique_GO_list(gene_inside_qtl_dict)
		
		#Prepare for counted genes inside qtl:
		golist_in_flat = flatten_array(gene_inside_qtl_dict)
		c_go_in_qtl_dict = count_all_goterms(golist_in_flat)
		
		#Prepare for counted genes outside qtl:
		golist_out_flat = flatten_array(gene_outside_qtl_dict)
		c_go_out_qtl_dict = count_all_goterms(golist_out_flat)
		
		#Create contingency tables for the Fishers exact test						
		c_array, total_genes = populate_contingency_table(c_go_in_qtl_dict, c_go_out_qtl_dict, golist_unique, tot_in_qtl, tot_out_qtl)
		
		#perform the fisher exact test on all created contingency tables
		fisher_python = fish_for_python(c_array)
		
		#only allow the results with p values below fisher_alpha to pass
		#fisher_alpha default is 0.05 unless another value is given
		#front end
		significant_info_fish = extract_significant_result_fish(fisher_python, fisher_alpha_name)
		
		#perform a multiple test
		enriched_golist_mult = correct_pvalues_for_multiple_testing(significant_info_fish)
		
		significant_info_mult = extract_significant_result_mult(enriched_golist_mult, mult_alpha_name)
		
		#Post processing
		qtl_go_dict = make_qtl_go_dict(gene_inside_qtl_dict, gene_dict)

		
		#A
		approved_golistA = post_process_A(enriched_golist_mult, qtl_go_dict, postA_name)
		#B
		godict_full = post_process_B(approved_golistA, gene_inside_qtl_dict, gene_outside_qtl_dict, postB_name)
		#go_gene_dict_full:
		#dict[i, go, fu_p_value, adj(fu_p_val, go_frac_scA, go_frac_scB]) = [gene list]

		
		del gene_list
		del absent_genes
		del present_genes
		del gene_inside_qtl_dict
		del gene_outside_qtl_dict
		del tot_in_qtl
		del tot_out_qtl
		del golist_unique
		del golist_in_flat
		del c_go_in_qtl_dict
		del golist_out_flat
		del c_go_out_qtl_dict
		del c_array
		del total_genes
		del fisher_python
		del significant_info_fish
		del significant_info_mult
		del enriched_golist_mult
		del qtl_go_dict
		del approved_golistA 

		return godict_full, qtls
	else:	
		return {}



		

def create_genelist_to_file(fn, dataset, cutoff):
	"""
	"""
	print "creating file for %s"%fn
	with open(fn, 'w') as fo:
		fo.write("dataset: %s"%dataset)
		fo.write("\n")
		fo.write("cutoff: %s"%cutoff)
		fo.write("\n\n")
		

def append_genelist_to_file(fn, trait, eqtl, genes):
	"""
	
	"""
	print "winner winner dicken chinner"
	print "writing data for %s"%trait
	with open(fn, 'a') as fo:
		fo.write("trait: %s "%trait)
		fo.write("\n")
		fo.write("eqtl: %s "%eqtl)
		fo.write("\n")
		for g in genes:
			fo.write("%s "%g)
		fo.write("\n\n")
		

def get_genes_from_godict(godict):
	"""
	"""
	genes = []
	for info, genelist in godict.iteritems():
		
		genes.extend(genelist)
	
	genes = sorted(list(set(genes)))
	return genes


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


def get_annotation_data():
	"""
	dont do this for each trait iteration...for God sakes...
	once is really enough!
	"""
	#get the location of the GO annotation file
	filename_annotate = 'At_geneGOlists.csv'
	module_dir = os.path.dirname('__file__')  # get current directory
	file_path_annotate = os.path.join(module_dir, 'qtl/documents', filename_annotate)
	
	#read the rows from the csv file mentioned above in filename
	data_dict = read_once_csv(file_path_annotate)
	
	return data_dict

	

def main():
	"""
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/get_enrichment.py')		
	"""

	#Datasets
	#exp_list = ['Snoek_2012']
	#exp_list = ['Ligterink_2014']
	exp_list = ['Keurentjes_2007']
	#Variables
	chromosome = [1,2,3,4,5]
	cutoff_list = [3]
	gxe_bool = False
	fisher_alpha_name = 0.05
	mult_alpha_name = 0.05
	postA_name = 0.5
	postB_name = 0.01
	
	anno_data = get_annotation_data()
	
	restart_necessary = True
	
	for dataset in exp_list:
		
		#make a new folder for the dataset
		if not gxe_bool:
			storage_folder = "%s/%s/enriched_%s"%(mr_folder, enriched_folder, dataset)
		else:
			storage_folder = "%s/%s/enriched_%s_gxe"%(mr_folder, enriched_folder, dataset)
		if not os.path.exists(storage_folder):
			os.mkdir(storage_folder)
		
		for cutoff in cutoff_list:
			
			traitlist = get_traits(dataset, cutoff, chromosome)
			##################################
			if restart_necessary:
				last_trait = "AT2G46230"
				trait_index = traitlist.index(last_trait)
				new_traitlist = traitlist[trait_index:]
			##################################
			print "--------------------------------------------------"
				
			if not gxe_bool:
				fname = "%s/enriched_%s_co%s.txt"%(storage_folder, dataset, cutoff)
			else:
				fname = "%s/enriched_%s_gxe_co%s.txt"%(storage_folder, dataset, cutoff)
			#print fname
			#prep file
			if not os.path.exists(fname):
				create_genelist_to_file(fname, dataset, cutoff)
				print "Created %s"%fname
				print "--------------------------------------------------"
			else:
				print "%s already exists"%fname
				print "--------------------------------------------------"
			if restart_necessary:
				
				for trait in new_traitlist:
					print trait	
					godict, eqtl = give_enrichment_for_trait(
											trait, cutoff, dataset, gxe_bool, 
											fisher_alpha_name, mult_alpha_name, 
											postA_name, postB_name, anno_data
										)
										
					genes = get_genes_from_godict(godict)
					
					append_genelist_to_file(fname, trait, eqtl, genes)
					
					del eqtl
					del genes
					gc.collect()
					
				del traitlist	
				del new_traitlist			
			else:
							
				for trait in traitlist:
					print trait	
					godict, eqtl = give_enrichment_for_trait(
											trait, cutoff, dataset, gxe_bool, 
											fisher_alpha_name, mult_alpha_name, 
											postA_name, postB_name, anno_data
										)

					genes = get_genes_from_godict(godict)
					
					append_genelist_to_file(fname, trait, eqtl, genes)
			
					del godict
					del eqtl
					del genes	
					gc.collect()
					
				del traitlist
main()




					
