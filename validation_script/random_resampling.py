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

import folder_assignments as fa
from data_handlers import (
			read_data, parse_AtReg_data, TFloc_pairs_AtRegNet, 
			parse_verdict
			)
			
from confusion_matrix import (
			identify_true_false_positives, count_false_negatives, 
			calculate_confusion, print_results
			)



def get_trait_samplesize_data(data):
	"""
	Data was written to:
	Home/xenv_dj1.6/QTL/validation/validation_numerics/tt_te_combi.txt
	
	This is data is to be extracted and will be partitioned based on
	dataset and cutoff
	
	"""

	data_dict = OrderedDict()
	
	for line in data:
		info = line.split()
		dataset = info[0].strip()
		cutoff = float(info[1].strip())
		trait = info[2].strip()
		sample_size = int(info[3].strip())
		
		if (dataset, cutoff) not in data_dict:
			data_dict[(dataset, cutoff)] = [[trait, sample_size]]
		else:
			temp_data = data_dict[(dataset, cutoff)]
			temp_data.append([trait, sample_size])
			data_dict[(dataset, cutoff)] = temp_data
			
	return data_dict


def get_genelist(dataset, cutoff, chromosome):
	"""
	Extract data from several files
	The files are arranged per dataset and cutoff
	
	The data consists of traits with a genelist from their eQTL regions
	"""
	
	data = []
	
	subfolder = "genelist_%s/"%dataset
			
	#Per trait extract the eQTL genelist
	for chrom in chromosome:
	
		if dataset == 'Ligterink_2014_gxe':
			datset = 'Ligterink_2014'

			filelocation = "%s/%s/genelist_%s/genelist_%s_co%s_chr%s.txt"%(
																	fa.mr_folder, 
																	fa.gfolder, datset, 
																	dataset, cutoff, chrom
																	)
		
		else:
			
			filelocation = "%s/%s/%s/genelist_%s_co%s_chr%s.txt"%(
																fa.mr_folder, 
																fa.gfolder, 
																subfolder, dataset, 
																cutoff, chrom
																)
																
		#These files contain the genelists from the eQTLs of every trait
		with open(filelocation, 'r') as fo:
			for line in fo:
				
				if line.startswith('trait'):
					trait = line[7:16]
					
				if line.startswith('AT'):
					genelist = line.split()

					if trait and genelist:
						data.append([trait, genelist])
				
	return data



def read_seeds(fn):
	seed_data = read_data(fn)
	return seed_data

	
	

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


def sum_of_TGTF_rel(trait_samplelist, tt_genes, reference):
	"""
	"""
	TGTF_rel_major = []
	
	for tr_ge in trait_samplelist:
		trait, samplelist = tr_ge
		g_tempset = set(samplelist)
		
		TF_inlist = g_tempset & reference	
	
		TGTF_rel_minor = identify_prediction(trait, samplelist, tt_genes, reference)
		TGTF_rel_major.extend(TGTF_rel_minor)
	
	return TGTF_rel_major


def identify_prediction(trait, genelist, true_target_genes, reference):
	"""
	Take a trait and a genelist
	And compare it to a reference regulatory network (AtRegNet)
	"""
	TG_TF_predictions = []
	
	e_tempset = set(genelist)
	#TFs that are found in the genelist
	TF_in_genelist = e_tempset & reference	
	#only copy the TG-TF relations for the traits that
	#were identified as having a TF in the genelist
	if TF_in_genelist and trait in true_target_genes:
		
		for e_id_TF in TF_in_genelist:
			
			#These are the verified TG-TF connections
			TG_TF_predictions.append([trait, e_id_TF])


	return TG_TF_predictions



def read_predicted_confusion_data(fn):
	"""
	"""
	data_dict = OrderedDict()
	
	conf_data = read_data(fn)
	
	for line in conf_data:
		if line.startswith("dataset:"):
			dataset = line[9:].strip()
		if line.startswith("cutoff:"):
			cutoff = float(line[8:].strip())
		if line.startswith("recall:"):
			recall = float(line[8:].strip())
		if line.startswith("precision:"):
			precision = float(line[11:].strip())
			
			if dataset and cutoff and recall and precision:
				data_dict[(dataset, cutoff)] = [recall, precision]
			
	return data_dict
	
	
			


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
	exp_list = ['Snoek_2012']#'Ligterink_2014','Ligterink_2014_gxe','Keurentjes_2007',
	cutoff_list = [6.7]#, 4.3, 6.7]
	chromo = [1,2,3,4,5]
	#-----------------------------------------------------------
	
	#compare every recall and precision that are calculated during the
	#random resampling to the ones in standard_data_dict per (dataset, cutoff)
	standardfileloc = "%s/%s/predicted_confusion_variables.txt"%(fa.mr_folder, fa.numfolder)
	standard_ref_data_dict = read_predicted_confusion_data(standardfileloc)
	
	
	fileloc = "%s/%s/tt_te_combi.txt"%(fa.mr_folder, fa.numfolder)
	szdata = read_data(fileloc)
	sample_size_dict = get_trait_samplesize_data(szdata)
	
	
	#get the premade 1000 distinct seeds of 8 digits each
	seedfile = "%s/%s/random_seeds.txt"%(fa.mr_folder, fa.numfolder)
	seeds = read_seeds(seedfile)
	
	resultsfolder = "%s/%s/permutate_result_test_F1.txt"%(fa.mr_folder, fa.numfolder)
	
	data_dict = {}
	

	for dataset in exp_list:
		
		for cutoff in cutoff_list:
			
			key = (dataset, cutoff)	
			
			#trait_genelist_list = [trait, genelist]
			trait_genelist_list = get_genelist(dataset, cutoff, chromo)
			
			#total_rel is summed over all traits in a (dataset, cutoff) combination
			total_rel = sum_of_total_rel(
									trait_genelist_list, 
									sh_TG_list_ref, 
									sh_TF_list_ref
									)
			
			if key in standard_ref_data_dict:
				ref_recall, ref_precision = standard_ref_data_dict[key]
				ref_F1 = 2 * (ref_precision * ref_recall) / (ref_precision + ref_recall)
				print "ref_F1: %s"%ref_F1
				
			#higher_recall = 0
			#lower_recall = 0
			#higher_precision = 0
			#lower_precision = 0
			
			higher_F1 = 0
			lower_F1 = 0


			if key in sample_size_dict:
				#sample_size_list = [trait, sample_size]
				sample_size_list = sample_size_dict[key]
				tt_genes = [item[0] for item in sample_size_list]
				
				#Go through the trait_genelist_list
				#Check for every passing trait if it is in trait_list
				#if so:
				#activate the randomizor!
				#by giving it the genelist, samplesize and a seed
			
			
			
			#permutate!
			for seedling in seeds:
				
				trait_randomsample = []

				#create [trait - sample gene list]
				for tr_ge in trait_genelist_list:
					g_trait, g_genelist = tr_ge
					
					if g_trait in tt_genes:
						trait_index = tt_genes.index(g_trait)
						s_trait, s_size = sample_size_list[trait_index]
						
						rsamp = select_random_sample(g_genelist, s_size, seedling)
						trait_randomsample.append([s_trait, rsamp])



				#TG_TF_pred is summed over all traits in a (dataset, cutoff) combination
				TG_TF_pred = sum_of_TGTF_rel(trait_randomsample, tt_genes, TF_set_ref)

	
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
				TP, FP, FN, TN, recall, specif, precision = calculate_confusion(
											total_rel, true_pred_rel, 
											false_pred_rel, unpredicted_rel
											)
				

				###########################################################
				#print "true_traits: %s"%len(set(tt_genes))
				#print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precision)
				###########################################################
				
				if precision + recall != 0:
					F1 = 2 * (precision * recall) / (precision + recall)
					print F1
				else:
					F1 = 0
					print "Division by zero prevented for %s co %s"%(dataset, cutoff)
				
				if F1 <= ref_F1:
					lower_F1 += 1
				if F1 > ref_F1:
					higher_F1 += 1

				
				#if recall < ref_recall:
					#lower_recall += 1
				#if recall > ref_recall:
					#higher_recall += 1
				#if precision < ref_precision:
					#lower_precision += 1
				#if precision > ref_precision:
					#higher_precision += 1
					
				###########################################################
			#toc = time.clock()
			#stopwatch = toc-tic
				
			with open(resultsfolder, 'a') as fo:
				
				fo.write("-------------------------")
				fo.write("\n")
				#fo.write("time: %s"%stopwatch)
				#fo.write("\n")
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
				#fo.write("recall:")
				#fo.write("\n")
				#fo.write("lower: %s"%lower_recall)
				#fo.write("\n")
				#fo.write("higher: %s"%higher_recall)
				#fo.write("\n")
				#fo.write("precision:")
				#fo.write("\n")
				#fo.write("lower: %s"%lower_precision)
				#fo.write("\n")
				#fo.write("higher: %s"%higher_precision)
				#fo.write("\n")
				fo.write("-------------------------")
				fo.write("\n")
				fo.write("\n")
				
				
			#print "-------------------------"
			#print "time: %s"%stopwatch
			#print "dataset: %s"%dataset
			#print "cutoff: %s"%cutoff
			#print "-------------------------"
			#print "recall:"
			#print "lower: %s"%lower_recall
			#print "higher: %s"%higher_recall
			#print "precision:"
			#print "lower: %s"%lower_precision
			#print "higher: %s"%higher_precision
			#print "-------------------------"
			print "-------------------------"
			print "dataset: %s"%dataset
			print "cutoff: %s"%cutoff
			print "-------------------------"
			print "lower_F1:\t%s"%lower_F1
			print "higher_F1:\t%s"%higher_F1	
			print "-------------------------"
					
		




main()








