"""
NB: This module makes use of the Django databases and must therefor be
used separate from the other validation modules

It will generate enrichment output for all traits in the database for a
particular dataset and cutoff

To run this scrip do the following:
Enter this command to open the django shell 
python manage.py shell

And enter this command to run this script	
execfile('qtl/enrichment_results.py')
"""

import sys
import os
import gc

from sortedcontainers import SortedDict

#The following lines enable this script to contact django
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD,Experiment,ExperimentMarker

from qtl.genelist_from_eqtl import (
							marker_logp_list, add_chromosome_and_position, 
							retrieve_logp_above_cutoff, check_region, 
							iterate_marker_tuple, highest_tuple, 
							set_region_from_adjacent_markers, 
							get_chromosome_max, check_regions_on_chromosome, 
							combine_marker_region_data, find_genes_in_regions, 
							normalize_object_data, read_distinct_genes
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
from qtl.prepare_for_display import get_genes_with_go_from_qtl




def give_genelist_for_trait(trait, cutoff, exp, gxe_boolean, fisher_alpha_name, mult_alpha_name, postA_name, postB_name, data_dict):
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
		#del tot_in_qtl
		#del tot_out_qtl
		del golist_unique
		del golist_in_flat
		#del c_go_in_qtl_dict
		del golist_out_flat
		#del c_go_out_qtl_dict
		del c_array
		del total_genes
		del fisher_python
		del significant_info_fish
		del significant_info_mult
		del enriched_golist_mult
		del qtl_go_dict
		del approved_golistA 

		return godict_full, gene_dict, qtls, c_go_in_qtl_dict, c_go_out_qtl_dict, tot_in_qtl, tot_out_qtl
	
	else:	
		return {}, {}, qtls, {}, {}, 0, 0



def make_dict(gene_dict, go_gene_dict):
	"""
	Prepare the display of the gene list for each go term 
	Divide the genelist of each go term into presence in QTLs
	
	#What is in these dictionaries?
	gene_dict : dict[(trait, chr nmbr, marker start pos, marker end pos)] = [gene name, gene start, gene end] 
	go_gene_dict : dict[i, go, fu_p_value, adj(fu_p_val, go_frac_scA, go_frac_scB]) = [gene list]
	go_definition_dict : dict[go term] = [name, definition]
	
	#Now what is really in these dictionaries?
	gene_dict contains eQTLs and a list of genes for each eQTL
	go_gene_dict contains GO terms (with statistical results) and a list of genes that were enriched for said GO term
	
	NB go_gene_dict == godict_full
	
	#By merging gene_dict and go_gene_dict it is possible to arrange GO terms and enriched genes into eQTLS
	go_genes_qtl_dict : dict[(url link, GO term, name, definition)] = list of [link, gene, description]
	
	
	#Edit: This function has been copied from "prepare_for_display.py" and altered
	"""
	
	go_dict = SortedDict()

	
	for go_info in go_gene_dict:
		
		go = go_info[1]
		pval = go_info[2]
		qval = go_info[3]
		qtl_genelist_array = []			
		intersecting_array = []
		
		go_gene_set = set(go_gene_dict[go_info])
	
		for qtl in gene_dict:
			

			
			#Create string that indicates the location of an eQTL
			#by displaying chromosome number and physical positions
			eQTL_start = "%s_%s"%(qtl[1], qtl[2])
			eQTL_finish = "%s_%s"%(qtl[1], qtl[3])
			region = "%s-%s"%(eQTL_start, eQTL_finish)
			
			#...
			gene_info = gene_dict[qtl]
			gene_list = [i[0] for i in gene_info]
			qtl_gene_set = set(gene_list)
			
			#Another set() operation
			#the intersecting genes are the ones from an eQTL with a 
			#particular GO term
			intersecting_genes = go_gene_set & qtl_gene_set
			
			sorted_intersecting_genes = sorted(intersecting_genes)
			temptup = (region, pval, qval, sorted_intersecting_genes)
			qtl_genelist_array.append(temptup)
			
			
			del eQTL_start
			del eQTL_finish
			del region
			del gene_info
			del gene_list
			del qtl_gene_set
			del intersecting_genes
			del sorted_intersecting_genes
			del temptup


		go_dict[go] = qtl_genelist_array
		
		del go
		del qtl_genelist_array
		del intersecting_array	

	del gene_dict
	del go_gene_dict

	 
			
	#now we have a dictionary containing for each go term a list of genes
	#each list of genes is separated into occurence in respective eQTLs
	
	return go_dict



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



########################################################################
########################################################################
########################################################################



def prepare_file(fn, chrom, exp, cutoff):
	"""
	
	"""
	with open(fn, 'w') as fo:
		fo.write("chromosome %s"%chrom)
		fo.write("\n")
		fo.write("dataset %s"%exp)
		fo.write("\n")
		fo.write("cutoff %s"%cutoff)
		fo.write("\n\n")
		
		
def write_data(fn, trait, godict, cgoes_in, cgoes_out, eqtls, tot_in, tot_out):
	"""
	
	"""
	
	with open(fn, 'a') as fo:
		
		fo.write("trait: %s"%trait)
		fo.write("\n")
		fo.write("nr of eQTLS: %s"%eqtls)
		fo.write("\n")
		fo.write("tot genes in: %s"%tot_in)
		fo.write("\n")
		fo.write("tot genes out: %s"%tot_out)
		fo.write("\n")
					
		for go, eqtl in godict.iteritems():
			fo.write(go)
			fo.write("\n")
			
			for region, p, q, genelist in eqtl:
				fo.write(region)
				fo.write("\n")
				fo.write("pval: %s"%p)
				fo.write("\n")
				fo.write("qval: %s"%q)
				fo.write("\n")
				if go in cgoes_in:
					fo.write("go in: %s"%cgoes_in[go])
					fo.write("\n")
				
				if go in cgoes_out:
					fo.write("go out: %s"%cgoes_out[go])
					fo.write("\n")					
				
				for gene in genelist:
					fo.write("%s "%gene)
				fo.write("\n")
					
		fo.write("\n\n")
				
	



def print_dict(godict):
	"""
	testing procedures
	"""
	
	for go, eqtl in go_dict.iteritems():
		print go
		
		for region, genelist in eqtl:
			print region
			
			for gene in genelist:
				print gene

########################################################################
########################################################################



def enriched_results_for_all():
	"""
	ITerate over all traits in all chromosomes
	And store for each trait the genes that were enriched
	
	
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/enrichment_results.py')
	"""
	



	#Datasets
	exp_list = ['Ligterink_2014','Keurentjes_2007','Snoek_2012']

	#Variables
	chromosome = [1,2,3,4,5]
	cutoff_list = [6.7]#, 4.3, 3, 2]
	gxe_boolean = False
	fisher_alpha_name = 0.05
	mult_alpha_name = 0.05
	postA_name = 0.5
	postB_name = 0.01
	
	folder = "enrichment_results"
	
	anno_data = get_annotation_data()
	
	if not os.path.exists(folder):
		os.mkdir(folder)
	else:
		pass
		
	for exp in exp_list:
		if not os.path.exists(folder):
			os.mkdir(folder)

	##run Forest! run!
	
	for cutoff in cutoff_list:
		
		for dataset in exp_list:
			
			for chrom in chromosome:
				
				subfolder = "v4_%s_%s"%(dataset, cutoff)
				storage = "%s/%s"%(folder, subfolder)
				
				if not os.path.exists(storage):
					os.mkdir(storage)
				
				filename = "%s/%s/results4_%s_co%s_chr%s.txt"%(folder, subfolder, dataset, cutoff, chrom)
				
				if not os.path.exists(filename):
					prepare_file(filename, chrom, dataset, cutoff)
					print "Created %s"%filename
				else:
					print "%s already exists"%filename
			

				
				##################################
				traitfile = '/home/wouter/xenv_dj1.6/QTL/validation/trait_lists/traitlist_chr%s.txt'%chrom
				with open(traitfile, 'r') as fo:
					data = fo.readlines()
				
				traitlist = [trait.strip() for trait in data]
				
				#If the iteration ended, copypaste the last gene into LOCUS_ID
				#Change the name in the for loop to n_traitlist
				#The index will start from that point
				
				index = traitlist.index('AT1G67480')
				n_traitlist = traitlist[index:]
				
				for trait in n_traitlist:
				##################################
				
					go_gene_dict, gene_dict, qtls, counted_goes_in, counted_goes_out, tot_in, tot_out = give_genelist_for_trait(
										trait, cutoff, dataset, gxe_boolean, fisher_alpha_name, 
										mult_alpha_name, postA_name, postB_name, anno_data
										)
					
					if qtls != 0:
						print trait
						
						gdict = make_dict(gene_dict, go_gene_dict)
						write_data(filename, trait, gdict, counted_goes_in, counted_goes_out, qtls, tot_in, tot_out)
						
						del trait
						del gdict

					else:
						pass
			
					#print_dict(gdict)
					
		
					del go_gene_dict
					del gene_dict
					del qtls
					del counted_goes_in
					del counted_goes_out
					del tot_in
					del tot_out
								
					gc.collect()
					
					
				del data
				del traitlist

#make it so*
enriched_results_for_all()























#*make hand movement not unlike Jean Luc Picard while starting this script


