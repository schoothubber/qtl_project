"""
Take a random sample from the eQTL genelist for each trait

the random sample size is equal to the actual corresponding enrichment size

for these random samples perform the steps to create confusion matrices

calculate the F1 score with the recalls and precisions

compare the F1 scores
"""
from collections import OrderedDict
import random
import time
import os

import folder_assignments as fa
from data_handlers import (
			read_data, parse_AtReg_data, TFloc_pairs_AtRegNet, 
			parse_verdict, read_seeds, get_trait_samplesize_data
			)
			
from confusion_matrix import (
			identify_true_false_positives, count_false_negatives, 
			calculate_confusion, print_results
			)

from process_datapoints import (
			get_TGTF_from_genelist, count_total_relations, process_enrichment
			)


def select_random_sample(genelist, sample_size, seedling):
	"""
	Take a list of genes
	shuffle them
	and return a random sample of genes
	
	The random genes are selected from the first occurring genes in
	a shuffled genelist
	"""
	random.seed(seedling)
	random.shuffle(genelist)
	random_sample = genelist[0:sample_size]
	
	return random_sample
	

def sum_of_total_rel(trait_genelist, TG_list_ref, TF_list_ref):
	"""
	This function calculates the total number of TG-TF connections 
	present in the full eQTL genelist, per trait, per dataset-cutoff
	"""
	total_rel_major = []
	refset = set(TF_list_ref)
	
	for tr_ge in trait_genelist:
		trait, genelist = tr_ge
		g_tempset = set(genelist)
		
		TF_in_genelist = g_tempset & refset
		
		if TF_in_genelist:

			total_rel_minor = count_total_relations(
											trait, TF_in_genelist, 
											TG_list_ref, TF_list_ref
											)
			total_rel_major.extend(total_rel_minor)
		
		
	
	return total_rel_major


def count_total_relations(trait, TF_in_genelist, TG_list_ref, TF_list_ref):
	"""
	
	"""
	total_relations = []
	if TF_in_genelist:
		for g_id_TF in TF_in_genelist:				
				
			if trait in TG_list_ref and g_id_TF in TF_list_ref:
				total_relations.append([trait, g_id_TF])
				
	return total_relations



def get_randomized_predictions(trait_randomsample, TF_set_ref):
	"""
	"""
	TGTF_list = []
	for info in trait_randomsample:
		TG = info[0]
		genelist = info[1]
		geneset = set(genelist)
		
		TF_in_sample = geneset & TF_set_ref
		if TF_in_sample:
			for TF in TF_in_sample:
				TGTF_list.append([TG, TF])
	return TGTF_list



def read_predicted_confusion_data(fn):
	"""
	"""
	
	conf_data = read_data(fn)
	
	for line in conf_data:
		if line.startswith("dataset:"):
			dataset = line[9:].strip()
		if line.startswith("cutoff:"):
			cutoff = float(line[8:].strip())
			
		if line.startswith("recall"):
			recall = line[7:].strip()
			if recall != "None":
				recall = float(recall)
			else:
				recall = None
				
		if line.startswith("precis"):
			precision = line[7:].strip()
			if precision != "None":
				precision = float(precision)
			else:
				precision = None
		
		if line.startswith("F1"):
			F1 = line[3:].strip()
			if F1 != "None":
				F1 = float(F1)
			else:
				F1 = None
				
			if F1:
				return recall, precision, F1
	
	
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
	datalist_eqtl = []
	datadict = {}
	data = read_data(fn)
	for line in data:
		if line.startswith("trait:"):
			trait = line[7:].strip()
		if line.startswith("eqtl:"):
			eqtl = int(line[6:].strip())
		if line.startswith("AT"):
			genelist = line.split()
			
			datalist_eqtl.append([trait, eqtl])
			datadict[trait] = genelist			
			
	return datalist_eqtl, datadict


def get_truetraits(fn):
	"""
	"""
	truetrait_list = []
	data = read_data(fn)
	for line in data:	
		if line.startswith("AT"):
			trait = line.strip()
			truetrait_list.append(trait)
			
	return truetrait_list


def main():
	"""
	"""
	tic = time.clock()
	#Get test data from AtRegNet.txt
	AtReg_data = read_data(fa.filename_atreg)
	AtRegNet_parse = parse_AtReg_data(AtReg_data)
	
	TF_TG_ref = TFloc_pairs_AtRegNet(AtRegNet_parse, ['all'])
	TG_TF_ref = [[info[1], info[0]] for info in TF_TG_ref]
	
	TG_list_ref = [info[0] for info in TG_TF_ref]
	TG_set_ref = set(TG_list_ref)
	sh_TG_list_ref = list(TG_set_ref)
	
	TF_list_ref = [info[1] for info in TG_TF_ref]
	TF_set_ref = set(TF_list_ref)
	sh_TF_list_ref = list(TF_set_ref)

	#-----------------------------------------------------------
	#exp_list = ['Snoek_2012']
	exp_list = ['Ligterink_2014','Ligterink_2014_gxe','Snoek_2012','Keurentjes_2007']
	cutoff_list = [3,4.3,6.7]
	chromo = [1,2,3,4,5]
	#-----------------------------------------------------------
	
	eQTL_threshold_list = [0,1,2,3]
	
	#-----------------------------------------------------------
	
	#get the premade 1000 distinct seeds of 8 digits each
	seedfile = "%s/%s/random_seeds.txt"%(fa.mr_folder, fa.numfolder)
	#print "Retrieving randomized seeds from %s"%seedfile
	seeds = read_seeds(seedfile)
	
	data_dict = {}
	
	write_summary = False
	write_conf = False
	print_conf = True
	

	summary = []
	print "start:"
	print "----------------------------"

	for eQTL_threshold in eQTL_threshold_list:

		for dataset in exp_list:
			
			for cutoff in cutoff_list:
			

				#print "Initializing analysis for dataset %s with cutoff %s"%(dataset, cutoff)
				
				F1 = None
				ref_F1 = None
							
				############################################################		
				####Retrieve original confusion matrix results
				subfolder_F1 = "/eqtl_%s/valnum_%s"%(eQTL_threshold, dataset)
				F1_fn = "%s/%s/%s/valnum_results_%s_co%s"%(
												fa.mr_folder, fa.numfolder,
												subfolder_F1, dataset, cutoff
												)
				try:
					ref_recall, ref_precision, ref_F1 = read_predicted_confusion_data(F1_fn)
				except:
					ref_recall= ref_precision= ref_F1 = None
				
				if ref_F1 != None:
		
					############################################################
					####Retrieve genelist
					subfolder_genelist = "genelist_%s"%dataset
					genelist_fn = "%s/%s/%s/genelist_%s_co%s.txt"%(
													fa.mr_folder, fa.gfolder,
													subfolder_genelist, dataset,
													cutoff
													)
					trait_genelist_list = get_genelist(genelist_fn)
					true_rel, total_rel = get_TGTF_from_genelist(
															genelist_fn, TG_TF_ref, 
															TG_list_ref, TF_list_ref
															)
					tt_genes = list(set([info[0] for info in true_rel]))
					############################################################			
					####Retrieve enriched list
					subfolder_enriched = "enriched_%s"%dataset
					enriched_fn = "%s/%s/%s/enriched_%s_co%s.txt"%(
													fa.mr_folder, fa.enriched_folder,
													subfolder_enriched, dataset,
													cutoff
													)
					trait_eqtl_genelist, dict_trait_enriched = get_enriched(enriched_fn)
					
					############################################################
					####Retrieve True Traits
					emr_traits_fn = "%s/%s/emr_traitlist_%s_co%s.txt"%(
													fa.mr_folder, fa.trait_folder,
													dataset, cutoff
													)

		
					############################################################
					#get all traits that have more than X eQTLs, where X = eQTL_threshold			
					truetrait_eqtl_list = [[t[0], t[1]] for t in trait_eqtl_genelist if t[0] in tt_genes and t[1]>eQTL_threshold]
					trait_with_eqtl = [info[0] for info in truetrait_eqtl_list]
					

					print "dataset:", dataset
					print "cutoff:", cutoff
					print "eQTL_threshold:", eQTL_threshold
					print "traits with eQTL:", len(trait_with_eqtl)





main()





					
