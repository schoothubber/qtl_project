"""
The construction of artificial datasets is done so they can be used to
verify that the confusion matrix construction is correctly executed
"""
import copy
import random

import folder_assignments as fa
from data_handlers import (
			read_data, parse_AtReg_data, TFloc_pairs_AtRegNet, parse_verdict
			)
from confusion_matrix import (
			identify_true_false_positives, count_false_negatives, 
			calculate_confusion, print_results
			)




def add_false_positives(TGTF_ad, TG_list, TF_list, perc):
	"""
	False Positives are those blue dots that are predicted TFs but
	are not verified by the reference data
	
	calculate howmany TFs should be added to the wrong TG
	"""
	dataset_size = len(TGTF_ad)
	temp_data = copy.deepcopy(TGTF_ad)
	
	if perc != 0:

		add_noise = perc * dataset_size / (100-perc)
		#
		#pick random TFs from AtRegNet and add them to a wrong TG
		#This part adds False Positives; Blue dots
		#
		#Each percentage should reflected by the precision
		enough = False
		while not enough:
			rand_TG = random.choice(TG_list)
			rand_TF = random.choice(TF_list)

			
			if [rand_TG, rand_TF] not in TGTF_ad:
				
				temp_data.append([rand_TG, rand_TF])
				add_noise -= 1
			
			if add_noise == 0:
				enough = True
	else:
		return temp_data

	return temp_data
					



def add_false_negatives(TGTF_ad, TG_TF_ref, perc):
	"""
	False Negatives are the red dots that were not overlapped by blue dots
	Also their trait must be a true trait
	
	i.e. remove a number of true [TG-TF]s from the artificial dataset
	this becomes a reduced artificial dataset
	"""
	
	dataset_size = len(TGTF_ad)
	temp_data = copy.deepcopy(TGTF_ad)
	removal = dataset_size * (perc / float(100))
	
	while removal > 0:
		rand_TGTF = random.choice(TGTF_ad)
		if rand_TGTF in temp_data:
			temp_data.remove(rand_TGTF)
			removal -= 1

	
	return temp_data
				








def main():
	"""
	confusion_matrix.py:
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
	TP, FP, FN, TN, recall, specif, precis = calculate_confusion(
								total_rel, true_pred_rel, 
								false_pred_rel, unpredicted_rel
								)
	###########################################################
	print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precis)
	###########################################################
	"""
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
	###########################################################
	###########################################################
	artificial_data = copy.deepcopy(TG_TF_ref)
	dataset = "artificial dataset"
	cutoff = 0
	
	###########################################################
	###########################################################
	if true_P:
		tt_genes = [info[0] for info in artificial_data]
		total_rel = artificial_data
		###########################################################								
		true_pred_rel, false_pred_rel = identify_true_false_positives(
															artificial_data,
															TG_TF_ref
															)
		###########################################################

		unpredicted_rel = count_false_negatives(
												TG_TF_ref, true_pred_rel, 
												tt_genes
												)										
		###########################################################
		TP, FP, FN, TN, recall, specif, precis = calculate_confusion(
									total_rel, true_pred_rel, 
									false_pred_rel, unpredicted_rel
									)
		###########################################################
		print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precis)
	###########################################################	
	###########################################################


	percentages = [20, 40, 50, 60]



	if not false_N and false_P:
		#Add some noise in the form of False Positives
		for perc in percentages:
			cutoff = perc
			noised_artificial_data = add_false_positives(
										artificial_data, sh_TG_list_ref, 
										sh_TF_list_ref, perc
										)
			###########################################################
			tt_genes = [info[0] for info in artificial_data]
			total_rel = noised_artificial_data
			###########################################################								
			true_pred_rel, false_pred_rel = identify_true_false_positives(
																noised_artificial_data,
																TG_TF_ref
																)
			###########################################################

			unpredicted_rel = count_false_negatives(
													TG_TF_ref, true_pred_rel, 
													tt_genes
													)										
			###########################################################
			TP, FP, FN, TN, recall, specif, precis = calculate_confusion(
										total_rel, true_pred_rel, 
										false_pred_rel, unpredicted_rel
										)
			###########################################################
			print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precis)
		###########################################################	
		###########################################################		


	if false_N and not false_P:
		#Add some noise in the form of False Negatives
		for perc in percentages:
			cutoff = perc
			noised_artificial_data = add_false_negatives(artificial_data, TG_TF_ref, perc)
			
			###########################################################
			tt_genes = [info[0] for info in artificial_data]
			total_rel = artificial_data
			###########################################################								
			true_pred_rel, false_pred_rel = identify_true_false_positives(
																noised_artificial_data,
																TG_TF_ref
																)
			###########################################################

			unpredicted_rel = count_false_negatives(
													TG_TF_ref, true_pred_rel, 
													tt_genes
													)										
			###########################################################
			TP, FP, FN, TN, recall, specif, precis = calculate_confusion(
										total_rel, true_pred_rel, 
										false_pred_rel, unpredicted_rel
										)
			###########################################################
			print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precis)
		###########################################################	
		###########################################################	



	if false_P and false_N:
		#Add some noise in the form of False Positives
		#And then add some False_Negatives
		for perc in percentages:
			cutoff = perc
			print "ad:",len(artificial_data)
			noised_artificial_data = add_false_positives(
													artificial_data, 
													sh_TG_list_ref, 
													sh_TF_list_ref, perc
													)	
			print "nad:",len(noised_artificial_data)
			more_noised_artificial_data = add_false_negatives(
												noised_artificial_data, 
												TG_TF_ref, perc
												)
			print "mnad:",len(more_noised_artificial_data)
		
			
			
			###########################################################
			tt_genes = [info[0] for info in noised_artificial_data]
			total_rel = noised_artificial_data
			###########################################################								
			true_pred_rel, false_pred_rel = identify_true_false_positives(
																more_noised_artificial_data,
																TG_TF_ref
																)
			###########################################################

			unpredicted_rel = count_false_negatives(
													TG_TF_ref, true_pred_rel, 
													tt_genes
													)										
			###########################################################
			TP, FP, FN, TN, recall, specif, precis = calculate_confusion(
										total_rel, true_pred_rel, 
										false_pred_rel, unpredicted_rel
										)
			###########################################################
			print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precis)
		###########################################################	
		###########################################################			

		
		
#What to do?
true_P = False
false_N = True
false_P = True
		
		
		
main()




















