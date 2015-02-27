import folder_assignments as fa
from data_handlers import (
			read_data, parse_AtReg_data, TFloc_pairs_AtRegNet, parse_verdict
			)
			
from confusion_matrix import (
			identify_true_false_positives, count_false_negatives, 
			calculate_confusion, print_results
			)
			
from process_datapoints import (
			get_TGTF_from_genelist, count_total_relations, process_enrichment,
			identify_prediction
			)



def main():
	"""
	"""
	#Get test data from AtRegNet.txt
	AtReg_data = read_data(fa.filename_atreg)
	AtRegNet_parse = parse_AtReg_data(AtReg_data)
	
	TF_TG_ref = TFloc_pairs_AtRegNet(AtRegNet_parse, ['all'])
	TG_TF_ref = [[info[1], info[0]] for info in TF_TG_ref]
	
	TG_list_ref = [info[0] for info in TG_TF_ref]
	TG_set_ref = set(TG_list_ref)
	
	TF_list_ref = [info[1] for info in TG_TF_ref]
	TF_set_ref = set(TF_list_ref)

	#-----------------------------------------------------------
	exp_list = ['Ligterink_2014']#,'Ligterink_2014_gxe','Keurentjes_2007','Snoek_2012']
	cutoff_list = [3, 4.3, 6.7]
	chromo = [1,2,3,4,5]
	#-----------------------------------------------------------

	for dataset in exp_list:
		
		for cutoff in cutoff_list:
			
			print "Analysing %s %s"%(dataset, cutoff)
			
			###########################################################
			true_rel, total_rel = get_TGTF_from_genelist(
													dataset, cutoff, chromo, 
													TG_TF_ref, TG_list_ref, TF_list_ref
													)
			
			###########################################################
			tt_genes = list(set([info[0] for info in true_rel]))
			
			###########################################################
			TG_TF_pred = process_enrichment(
												dataset, cutoff, chromo, 
												TF_set_ref, tt_genes
											)
											
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
			print "true_traits: %s"%len(set(tt_genes))
			print_results(dataset, cutoff, TP, FP, FN, TN, recall, specif, precis)
			###########################################################					
			
main()
			





