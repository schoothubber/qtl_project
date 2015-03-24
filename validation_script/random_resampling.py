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
	exp_list = ['Snoek_2012']
	#exp_list = ['Ligterink_2014','Ligterink_2014_gxe','Snoek_2012','Keurentjes_2007']
	cutoff_list = [3]#[4.3,6.7]
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
	
	for dataset in exp_list:
		
		for cutoff in cutoff_list:
			
			for eQTL_threshold in eQTL_threshold_list:
				
				print "Initializing analysis for dataset %s with cutoff %s"%(dataset, cutoff)
				
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
					print "traits with eQTL", len(trait_with_eqtl)
					############################################################
					
					higher_recall=lower_recall=higher_precision=lower_precision=0
					higher_F1=lower_F1=0		
					
					permutated_confusion = []
					#permutate!
					#print "Commencing permutation of %s, standby..."%len(seeds)
					#i = 0
					for seedling in seeds:
						#reset variables
						TP=FP=FN=TN=recall=specif=precision=F1= 0
						#print i
						#i += 1
						
						trait_randomsample = []
		
						#create [trait - sample gene list]
						for tr_ge in trait_genelist_list:
							g_trait, g_genelist = tr_ge
							#print g_trait
							#print len(g_genelist)
							
							if g_trait in trait_with_eqtl and g_trait in dict_trait_enriched:
								
								sample_size = len(dict_trait_enriched[g_trait])
								rsamp = select_random_sample(g_genelist, sample_size, seedling)
								trait_randomsample.append([g_trait,0, rsamp])
								#q = len(g_genelist)
								#print "take %s from %s"%(sample_size, q)


						#TG_TF_pred is summed over all traits in a (dataset, cutoff) combination
						#TG_TF_pred = get_randomized_predictions(trait_randomsample, TF_set_ref)
						TG_TF_pred = process_enrichment(trait_randomsample, TF_set_ref)
		
						
						#proceed with the random sample to the confusion matrix
						###########################################################										
						true_pred_rel, false_pred_rel = identify_true_false_positives(
																			TG_TF_pred,
																			TG_TF_ref
																			)
						
						###########################################################
						unpredicted_rel = count_false_negatives(
																TG_TF_ref, true_pred_rel, 
																tt_genes
																)
																
						###########################################################
						TP, FP, FN, TN, recall, specif, precision, F1 = calculate_confusion(
													total_rel, true_pred_rel, 
													false_pred_rel, unpredicted_rel
													)
						
						permutated_confusion.append([TP, FP, FN, TN, recall, specif, precision, F1])
						###########################################################
						#print "true_traits: %s"%len(set(tt_genes))
						#print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precision)
						###########################################################
						
		
						if recall != None:
							if recall < ref_recall:
								lower_recall += 1
							if recall >= ref_recall:
								higher_recall += 1
						#else:
							#print "recall is None"
							#pass
							
						if precision != None:
							if precision < ref_precision:
								lower_precision += 1
							if precision >= ref_precision:
								higher_precision += 1
						#else:
							#print "precision is None"
							#pass
							
						if ref_F1 != None and F1 != None:
								if F1 < ref_F1:
									lower_F1 += 1
								if F1 >= ref_F1:
									higher_F1 += 1
		
					summary.append([dataset, cutoff, lower_recall, higher_recall, lower_precision, higher_precision, lower_F1, higher_F1])
						###########################################################
					
					
				if write_conf:
					substorage = "%s/%s/%s"%(fa.mr_folder, fa.numfolder, dataset)
					if not os.path.exists(substorage):
						os.mkdir(substorage)
						
					resultsfolder_conf = "%s/permutate_eqtl_%s_%s_co%s.txt"%(
														 substorage, eQTL_threshold, 
														 dataset, cutoff
														)

					try:	
						print "Writing to file %s"%resultsfolder_conf
						with open(resultsfolder_conf, 'w') as fo:
							fo.write("-------------------------")
							fo.write("\n")
							fo.write("dataset: %s"%dataset)
							fo.write("\n")
							fo.write("cutoff: %s"%cutoff)
							fo.write("\n")
							fo.write("-------------------------")
							fo.write("\n")
							fo.write("lower_F1: \t%s"%lower_F1)
							fo.write("\n")
							fo.write("higher_F1: \t%s"%higher_F1)
							fo.write("\n")
							fo.write("-------------------------")
							fo.write("\n")
							fo.write("lower_recall: %s"%lower_recall)
							fo.write("\n")
							fo.write("higher_recall: %s"%higher_recall)
							fo.write("\n")
							fo.write("lower_precision: %s"%lower_precision)
							fo.write("\n")
							fo.write("higher_precision: %s"%higher_precision)
							fo.write("\n")
							fo.write("-------------------------")
							fo.write("\n")
							for [TP, FP, FN, TN, recall, specif, precision, F1] in permutated_confusion:
								fo.write("-------------------------\n")
								fo.write("TP\t%s\tFN\t%s"%(TP, FN))
								fo.write("\n")
								fo.write("FP\t%s\tTN\t%s"%(FP, TN))
								fo.write("\n")
								fo.write("-------------------------\n")
								fo.write("recall\t%s"%recall)
								fo.write("\n")
								fo.write("specificity\t%s"%specif)
								fo.write("\n")
								fo.write("precision\t%s"%precision)
								fo.write("\n")
								fo.write("F1\t%s"%F1)
								fo.write("\n")
								fo.write("-------------------------\n")
					except:
						pass
											
				if print_conf:
					try:
						print "-------------------------"
						print "TP\t%s\tFN\t%s"%(TP, FN)
						print "FP\t%s\tTN\t%s"%(FP, TN)
						print "-------------------------"
						print "dataset: %s"%dataset
						print "cutoff: %s"%cutoff
						print "eQTL: %s"%eQTL_threshold
						print "-------------------------"
						print "lower_F1:\t%s"%lower_F1
						print "higher_F1:\t%s"%higher_F1	
						print "-------------------------"
						print "lower_recall: %s"%lower_recall
						print "higher_recall: %s"%higher_recall
						print "lower_precision: %s"%lower_precision
						print "higher_precision: %s"%higher_precision
						print "-------------------------"
					except:
						pass


						
			if write_summary:
				summfolder_conf = "%s/%s/permutate_summary_eqtl_%s.txt"%(
												fa.mr_folder, fa.numfolder,
												eQTL_threshold
												)
				try:
					with open(summfolder_conf, 'w') as fo:
						for dataset, cutoff, lower_recall, higher_recall, lower_precision, higher_precision, lower_F1, higher_F1 in summary:
							fo.write("-------------------------")
							fo.write("\n")
							fo.write("dataset: %s"%dataset)
							fo.write("\n")
							fo.write("cutoff: %s"%cutoff)
							fo.write("\n")
							fo.write("-------------------------")
							fo.write("\n")
							fo.write("lower_F1: \t%s"%lower_F1)
							fo.write("\n")
							fo.write("higher_F1: \t%s"%higher_F1)
							fo.write("\n")
							fo.write("-------------------------")
							fo.write("\n")
							fo.write("recall:")
							fo.write("\n")
							fo.write("lower: %s"%lower_recall)
							fo.write("\n")
							fo.write("higher: %s"%higher_recall)
							fo.write("\n")
							fo.write("precision:")
							fo.write("\n")
							fo.write("lower: %s"%lower_precision)
							fo.write("\n")
							fo.write("higher: %s"%higher_precision)
							fo.write("\n")
							fo.write("-------------------------")
							fo.write("\n")
				except:
					pass




main()








