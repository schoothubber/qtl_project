import folder_assignments as fa
from data_handlers import (
			read_data, parse_AtReg_data, TFloc_pairs_AtRegNet, parse_verdict
			)
			
from confusion_matrix import (
			identify_true_false_positives, count_false_negatives, 
			calculate_confusion, print_results
			)
			
from process_datapoints import (
			get_TGTF_from_genelist, count_total_relations, process_enrichment
			)


def get_info(fn):
	"""
	"""
	data = read_data(fn)
	trait_eqtl_genelist = []
	
	for line in data:
		if line.startswith("trait:"):
			trait = line[7:16]
		if line.startswith("eqtl:"):
			eqtl = int(line[6:].strip())
		if line.startswith("AT"):
			genelist = line.split()
		
			trait_eqtl_genelist.append([trait, eqtl, genelist])
		
	return trait_eqtl_genelist



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
	#exp_list = ['Ligterink_2014','Ligterink_2014_gxe','Keurentjes_2007','Snoek_2012']
	exp_list = ['Ligterink_2014']
	cutoff_list = [3]#, 4.3, 6.7]
	chromo = [1,2,3,4,5]
	#-----------------------------------------------------------

	for dataset in exp_list:
		
		for cutoff in cutoff_list:
			
			print "Analysing %s %s"%(dataset, cutoff)
			
			###########################################################
			#Extract the true TG-TF relations and the total possible
			#relations from the stored datafiles
			filelocation = "%s/%s/genelist_%s/genelist_%s_co%s.txt"%(
															fa.mr_folder, 
															fa.gfolder, 
															dataset, dataset, 
															cutoff
															)
			true_rel, total_rel = get_TGTF_from_genelist(
													filelocation, TG_TF_ref, 
													TG_list_ref, TF_list_ref
													)
			
			###########################################################
			#The TG in the true TG-TF relations are the true_traits (tt)
			#in this case named tt_genes
			tt_genes = list(set([info[0] for info in true_rel]))
			
########################################################################
			#Get for each true_trait the number of eQTLs
			enriched_fn = "%s/%s/enriched_%s/enriched_%s_co%s.txt"%(
												fa.mr_folder, fa.enriched_folder, 
												dataset, dataset, cutoff
												)
			trait_eqtl_genelist = get_info(enriched_fn)
			#Select true traits based on number of eQTLs
			tt_trait_eqtl_genelist = [[t[0], t[1], t[2]] for t in trait_eqtl_genelist if t[0] in tt_genes and t[1]>0]
			trait_with_eqtl = [info[0] for info in tt_trait_eqtl_genelist]
########################################################################
			
			###########################################################
			TG_TF_pred = process_enrichment(tt_trait_eqtl_genelist, TF_set_ref)
											
			###########################################################										
			true_pred_rel, false_pred_rel = identify_true_false_positives(
																TG_TF_pred,
																TG_TF_ref
																)
			
			###########################################################
			unpredicted_rel = count_false_negatives(
													TG_TF_ref, true_pred_rel, 
													trait_with_eqtl
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
			





