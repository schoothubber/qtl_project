"""
This is the main script from which all validation should flow

Before this script you should run qtl/enrichment_results.py
in order to collect the enrichment data per trait

"""
import collections as col

import folder_assignments as fa
from data_handlers import (
			read_data, parse_AtReg_data, parse_family_data, 
			TFloc_pairs_AtRegNet, parse_verdict, make_AtReg_dict,
			write_pairwise_validation
			)
from verify_enrichment_results import verify
from vis_regulator import (
			get_all_loci, get_loci_from_dataset, create_reg_dict, 
			link_data_in_array, draw_plot, draw_combi_plot
			)
from process_datapoints import (
			get_TGTF_from_genelist, count_total_relations, process_enrichment,
			identify_prediction
			)
from confusion_matrix import (
			identify_true_false_positives, count_false_negatives, 
			calculate_confusion, print_results
			)			


#----------------------------------------------------------------------
#Set all variables
#----------------------------------------------------------------------
#exp_list = ['Ligterink_2014']
#exp_list = ['Ligterink_2014_gxe']
#exp_list = ['Keurentjes_2007']
exp_list = ['Snoek_2012']

chromosome = [1,2,3,4,5]
#cutoff_list = [3]
#cutoff_list = [4.3]
cutoff_list = [6.7]

dataset = exp_list[0]
cutoff = cutoff_list[0]		

#----------------------------------------------------------------------
#Set all data files and lists
#----------------------------------------------------------------------

#read data from files, parse and process

TF_family_data = read_data(fa.filename_fam)
TF_fam_data = parse_family_data(TF_family_data)
TF_fam_list = sorted([info[1] for info in TF_fam_data])
TF_fam_set = set(TF_fam_list)

AtReg_data = read_data(fa.filename_atreg)
AtRegNet_parse = parse_AtReg_data(AtReg_data)
AtRegNet_list = [info[2] for info in AtRegNet_parse]
AtRegNet_set = set(AtRegNet_list)
AR_dict = make_AtReg_dict(AtRegNet_parse)

fam_selection = ["all"]#use "all" for all regulators
AtRegNet_pairs = TFloc_pairs_AtRegNet(AtRegNet_parse, fam_selection)

TF_TG_ref = TFloc_pairs_AtRegNet(AtRegNet_parse, ['all'])
TG_TF_ref = [[info[1], info[0]] for info in TF_TG_ref]

TG_list_ref = [info[0] for info in TG_TF_ref]
TG_set_ref = set(TG_list_ref)
sh_TG_list_ref = list(TG_set_ref)

TF_list_ref = [info[1] for info in TG_TF_ref]
TF_set_ref = set(TF_list_ref)
sh_TF_list_ref = list(TF_set_ref)


#create the necessary dicts to create some scatter plots
loc_dict, chr_len_list = get_all_loci(fa.mr_folder, fa.trait_folder, len(chromosome))
reg_dict, TF_chr_len_list = create_reg_dict(list(AtRegNet_set))

#highly selective data sets
flowerlist = ["AT5G60910", "AT2G22540", "AT5G03840", "AT1G69120", "AT5G61850", "AT1G24260", "AT5G03790", "AT3G61250", "AT4G24540", "AT2G03710", "AT1G26310"]
flowerset = set(flowerlist)


#----------------------------------------------------------------------
#Reading and parsing the enrichment results
#----------------------------------------------------------------------

	#option 1 for the reference
	#option 2 for one of the models (just TF list from AtRegNet)
	#option 2.2 for one of the models (pairwise matrix from AtRegNet)
	#option 3 for a combination of reference and model
	#option 3.1 for visualisation of Genelist from validate_numerics.py
	#option 3.2 for visualisation of Enrichment from validate_numerics.py
	#option 4 for a combination of two models (Keurentjes vs Snoek)
	#option 5 for a more selective procedure
	#option 6 for a comparison between normal expression and gxe (Ligterink)

option = 3.1



if option == 1:
	
	
	fam_selection = ["all"]#use "all" for all regulators
	
	AtRegNet_pairs = TFloc_pairs_AtRegNet(AtRegNet_parse, fam_selection)
	print "AtRegNet_pairs %s"%len(AtRegNet_pairs)
	
	data_array, labels = link_data_in_array(AtRegNet_pairs, loc_dict, reg_dict, {})

	#Draw a scatter plot
	name = "AtRegNet"
	plot_title = "Gene regulation from %s"%name
	scatterfilelocation = "%s/%s/%s_plot3.png"%(fa.mr_folder, fa.scatplot_folder, name)
	print "------------------------------------"
	print "drawing %s"%scatterfilelocation
	print "------------------------------------"
	color = 'red'
	draw_plot(
				data_array, scatterfilelocation, plot_title, chr_len_list, 
				TF_chr_len_list, [], color
				)
	
	
	
	
if option == 2:
	
	for cutoff in cutoff_list:
			
		for dataset in exp_list:
			
			subfolder = "%s_%s/"%(dataset, cutoff)	
			combined_verdict = col.OrderedDict()
			print "------------------------------------"
			print "Analysing %s..."%subfolder
			print "------------------------------------"
			
			for chrom in chromosome:
				print "processing chromosome %s"%(chrom)
					
				filelocation = "%s/%s/v3_%s/results3_%s_co%s_chr%s.txt"%(
													fa.main_folder, 
													fa.enriched_folder, subfolder, 
													dataset, cutoff, 
													chrom
													)
				verdict = verify(AtRegNet_set, filelocation)
				#print "AtRegNet_set %s"%len(AtRegNet_set)
				
				#Combine dictionaries
				combined_verdict.update(verdict)
			
			print "------------------------------------"
			print "retrieving data for scatter plot..."
			print "------------------------------------"
		
			
			#From verdict extract the paired TF-loc list								
			TFloc_list, loc_list = parse_verdict(combined_verdict)
				
			#Prepare data for scatter plots
			data_array, label_list = link_data_in_array(TFloc_list, loc_dict, reg_dict, {})
			
	
			#Draw some scatter plots
			plot_title = "Predicted regulation from %s with cutoff %s"%(dataset, cutoff)
			scatterfilelocation = "%s/%s/scatter3_overlap_%s_co%s.png"%(
												fa.mr_folder, fa.scatplot_folder, 
												dataset, cutoff
												)
			print "------------------------------------"
			print "drawing %s"%scatterfilelocation
			print "------------------------------------"
			color = 'blue'
			draw_plot(
						data_array, scatterfilelocation, plot_title, 
						chr_len_list, TF_chr_len_list, [], color
						)
	
	
if option == 2.2:
	
	for cutoff in cutoff_list:
			
		for dataset in exp_list:
			
			subfolder = "%s_%s/"%(dataset, cutoff)	
			combined_verdict = col.OrderedDict()
			print "------------------------------------"
			print "Analysing %s..."%subfolder
			print "------------------------------------"
			
			for chrom in chromosome:
				print "processing chromosome %s"%(chrom)
				
				#add _gxe in front of v3_ to access the environmental perturbed dataset of Ligterink
				filelocation = "%s/%s/v3_%s/results3_%s_co%s_chr%s.txt"%(fa.main_folder, fa.enriched_folder, subfolder, dataset, cutoff, chrom)
				verdict = verify(AtRegNet_set, filelocation)
				#print "AtRegNet_set %s"%len(AtRegNet_set)
				
				#Combine dictionaries
				combined_verdict.update(verdict)
			
			print "------------------------------------"
			print "retrieving data for scatter plot..."
			print "------------------------------------"
		
			
			#From verdict extract the paired TF-loc list								
			TFloc_list, loc_list = parse_verdict(combined_verdict)
			
			#The loc_dict can be made from a any number of genes, these genes will form the x-axis in the scatterplot
			loc_dict, chr_len_list = get_loci_from_dataset(mr_folder, trait_folder, list(loc_list))
				
			#Prepare data for scatter plots
			data_array, label_list = link_data_in_array(TFloc_list, loc_dict, reg_dict, AR_dict)
			
			#Get the validated pairwise TF-gene relationships and write to a file
			val_pairwise_mod = [[row[2],row[3]] for row in label_list]
			
			pairwise_validationfilelocation = "%s/%s/network_%s_co%s"%(fa.mr_folder, fa.pvfolder, dataset, cutoff)
			write_pairwise_validation(pairwise_validationfilelocation, val_pairwise_mod)
			
			
	
			#Draw some scatter plots
			plot_title = "Predicted regulation from %s with cutoff %s"%(dataset, cutoff)
			scatterfilelocation = "%s/%s/scatter3_%s_co%s.png"%(fa.mr_folder, fa.scatplot_folder, dataset, cutoff)
			print "------------------------------------"
			print "drawing %s"%scatterfilelocation
			print "------------------------------------"
			draw_plot(data_array, scatterfilelocation, plot_title, chr_len_list, TF_chr_len_list, label_list)


			
if option == 3:
	
#Prepare data for scatter plots
	
	#get the reference data:
	
	print "------------------------------------"
	print "Retrieving data from AtRegNet.txt..."
	print "------------------------------------"
	
	AtRegNet_Data = read_data(filename_atreg)
	AtRegNet_parse = parse_AtReg_data(AtRegNet_Data)
	
	fam_selection = ["all"]#use "all" for all regulators
	AtRegNet_pairs = TFloc_pairs_AtRegNet(AtRegNet_parse, fam_selection)
	regdata_array, labels_reg = link_data_in_array(AtRegNet_pairs, loc_dict, reg_dict, {})

	
	#get the model data:
	for cutoff in cutoff_list:
			
		for dataset in exp_list:
			
			subfolder = "%s_%s"%(dataset, cutoff)	
			combined_verdict = col.OrderedDict()
			print "------------------------------------"
			print "Analysing %s..."%subfolder
			print "------------------------------------"
			
			for chrom in chromosome:
				print "processing chromosome %s"%(chrom)
					
				filelocation = "%s/%s/%s/results2_%s_co%s_chr%s.txt"%(fa.main_folder, fa.enriched_folder, subfolder, dataset, cutoff, chrom)
				verdict = verify(TF_fam_set, filelocation)
				#Combine dictionaries
				combined_verdict.update(verdict)
			
			print "------------------------------------"
			print "retrieving data for scatter plot..."
			print "------------------------------------"
		
				
			#From verdict extract the paired TF-loc list								
			TFloc_list, loc_list = parse_verdict(combined_verdict)
	
			#Prepare data for scatter plots
			moddata_array, labels_mod = link_data_in_array(TFloc_list, loc_dict, reg_dict, {})
	
			#Draw some scatter plots
			plot_title = "Predicted regulation from %s with cutoff %s"%(dataset, cutoff)
			scatterfilelocation = "%s/%s/scatter3_combined_%s_co%s.png"%(fa.mr_folder, fa.scatplot_folder, dataset, cutoff)
			print "------------------------------------"
			print "drawing %s"%scatterfilelocation
			print "------------------------------------"
			draw_combi_plot(moddata_array, regdata_array, scatterfilelocation, plot_title, chr_len_list, TF_chr_len_list)
	
	
			
if option == 3.1:
	#validate_validation() returns
	#(dataset, cutoff),temp_dict, TG_TF, truetraits, total
	print "------------------------------------"
	print "Retrieving data from AtRegNet.txt..."
	print "------------------------------------"
	
	label_lever = False
	color = 'blue'
	###########################################################		
	true_relations, total_rel_major = get_TGTF_from_genelist(
												dataset, cutoff, 
												chromosome, TG_TF_ref, 
												sh_TG_list_ref, sh_TF_list_ref
												)
	###########################################################	
	tt_genes = list(set([info[0] for info in true_relations]))
	###########################################################	
	TG_TF_pred = process_enrichment(
									dataset, cutoff, chromosome, 
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

	#TFloc_list contains predicted datapoints
	#TFloc_list = [[info[1],info[0]] for info in total_rel_major]
	
	regus_list = sorted(list(set([info[1] for info in total_rel_major])))
	locus_list = sorted(list(set([info[0] for info in total_rel_major])))

	regus_dict, TF_chr_len_list = create_reg_dict(regus_list)
	
	locus_dict, chr_len_list = get_loci_from_dataset(locus_list)
	
	#Prepare data for scatter plots
	moddata_array, labels_mod = link_data_in_array(total_rel_major, locus_dict, regus_dict, {})
	
	#Draw some scatter plots
	plot_title = "TG-TF in genelists of %s with cutoff %s"%(dataset, cutoff)
	scatterfilelocation = "%s/%s/scatter2_genelist_%s_co%s.png"%(
														fa.mr_folder, fa.scatplot_folder, 
														dataset, cutoff
														)
	print "------------------------------------"
	print "drawing %s"%scatterfilelocation
	print "------------------------------------"

	draw_plot(
		moddata_array, scatterfilelocation, plot_title, chr_len_list, 
		TF_chr_len_list, labels_mod, label_lever, color
		)
		
if option == 3.2:
	#validate_validation() returns
	#(dataset, cutoff),temp_dict, TG_TF, truetraits, total
	print "------------------------------------"
	print "Retrieving data from AtRegNet.txt..."
	print "------------------------------------"
	
	
	fam_selection = ["all"]#use "all" for all regulators
	AtRegNet_pairs = TFloc_pairs_AtRegNet(AtRegNet_parse, fam_selection)

	TF_TG_ref = TFloc_pairs_AtRegNet(AtRegNet_parse, ['all'])
	TG_TF_ref = [[info[1], info[0]] for info in TF_TG_ref]
	
	TG_list_ref = [info[0] for info in TG_TF_ref]
	TG_set_ref = set(TG_list_ref)
	sh_TG_list_ref = list(TG_set_ref)
	
	TF_list_ref = [info[1] for info in TG_TF_ref]
	TF_set_ref = set(TF_list_ref)
	sh_TF_list_ref = list(TF_set_ref)

	label_lever = False
	###########################################################		
	true_relations, total_rel_major = get_TGTF_from_genelist(
												dataset, cutoff, 
												chromosome, TG_TF_ref, 
												sh_TG_list_ref, sh_TF_list_ref
												)
	###########################################################	
	tt_genes = list(set([info[0] for info in true_relations]))
	###########################################################	
	TG_TF_pred = process_enrichment(
									dataset, cutoff, chromosome, 
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
	

	#TFloc_list contains predicted datapoints
	TFloc_list = [[info[1],info[0]] for info in true_pred_rel]
	
	regus_list = sorted(list(set([info[1] for info in total_rel_major])))
	locus_list = sorted(list(set([info[0] for info in total_rel_major])))

	regus_dict, TF_chr_len_list = create_reg_dict(regus_list)
	
	locus_dict, chr_len_list = get_loci_from_dataset(locus_list)
	
	#Prepare data for scatter plots
	moddata_array, labels_mod = link_data_in_array(TFloc_list, locus_dict, regus_dict, {})
	
	TF_TG_unrel = [[info[1], info[0]] for info in unpredicted_rel]
	regdata_array, labels_reg = link_data_in_array(TF_TG_unrel, locus_dict, regus_dict, {})
	
	#Draw some scatter plots
	plot_title = "TG-TF in Enrichment of %s with cutoff %s"%(dataset, cutoff)
	scatterfilelocation = "%s/%s/scatter_Enrichment_combined_%s_co%s.png"%(
														fa.mr_folder, fa.scatplot_folder, 
														dataset, cutoff
														)
	print "------------------------------------"
	print "drawing %s"%scatterfilelocation
	print "------------------------------------"
	draw_combi_plot(
				moddata_array, regdata_array, scatterfilelocation, 
				plot_title, chr_len_list, TF_chr_len_list, labels_mod,
				label_lever
				)
		


if option == 3.3:
	#validate_validation() returns
	#(dataset, cutoff),temp_dict, TG_TF, truetraits, total
	print "------------------------------------"
	print "Retrieving data from AtRegNet.txt..."
	print "------------------------------------"
	
	AtRegNet_Data = read_data(filename_atreg)
	AtRegNet_parse = parse_AtReg_data(AtRegNet_Data)
	
	fam_selection = ["all"]#use "all" for all regulators
	AtRegNet_pairs = TFloc_pairs_AtRegNet(AtRegNet_parse, fam_selection)

	TF_TG_ref = TFloc_pairs_AtRegNet(AtRegNet_parse, ['all'])
	TG_TF_ref = [[info[1], info[0]] for info in TF_TG_ref]
	
	TG_list_ref = [info[0] for info in TG_TF_ref]
	TG_set_ref = set(TG_list_ref)
	sh_TG_list_ref = list(TG_set_ref)
	
	TF_list_ref = [info[1] for info in TG_TF_ref]
	TF_set_ref = set(TF_list_ref)
	sh_TF_list_ref = list(TF_set_ref)

	label_lever = False
	color = 'blue'
	###########################################################		
	true_relations, total_rel_major = get_TGTF_from_genelist(
												dataset, cutoff, 
												chromosome, TG_TF_ref, 
												sh_TG_list_ref, sh_TF_list_ref
												)
	###########################################################	
	tt_genes = list(set([info[0] for info in true_relations]))
	###########################################################	
	TG_TF_pred = process_enrichment(
									dataset, cutoff, chromosome, 
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

	#TFloc_list contains predicted datapoints
	TFloc_list = [[info[1],info[0]] for info in total_rel_major]
	
	regus_list = sorted(list(set([info[1] for info in total_rel_major])))
	locus_list = sorted(list(set([info[0] for info in total_rel_major])))

	regus_dict, TF_chr_len_list = create_reg_dict(regus_list)
	
	locus_dict, chr_len_list = get_loci_from_dataset(locus_list)
	
	#Prepare data for scatter plots
	moddata_array, labels_mod = link_data_in_array(TFloc_list, locus_dict, regus_dict, {})
	
	#Draw some scatter plots
	plot_title = "TG-TF in genelists of %s with cutoff %s"%(dataset, cutoff)
	scatterfilelocation = "%s/%s/scatter_genelist_%s_co%s.png"%(
														fa.mr_folder, fa.scatplot_folder, 
														dataset, cutoff
														)
	print "------------------------------------"
	print "drawing %s"%scatterfilelocation
	print "------------------------------------"

	draw_plot(
		moddata_array, scatterfilelocation, plot_title, chr_len_list, 
		TF_chr_len_list, labels_mod, label_lever, color
		)




			
if option == 4:
	
	dataset_list = []
	#Get data from model 1 and 2
	for cutoff in cutoff_list:
			
		for dataset in exp_list:
			
			subfolder = "%s_%s"%(dataset, cutoff)	
			combined_verdict = col.OrderedDict()
			print "------------------------------------"
			print "Analysing %s..."%subfolder
			print "------------------------------------"
			
			for chrom in chromosome:
				print "processing chromosome %s"%(chrom)
					
				filelocation = "%s/%s/%s/results_%s_co%s_chr%s.txt"%(fa.mr_folder, fa.enriched_folder, subfolder, dataset, cutoff, chrom)
				verdict = verify(TF_fam_set, filelocation)
				#Combine dictionaries
				combined_verdict.update(verdict)
			
			print "------------------------------------"
			print "retrieving %s data for scatter plot..."%subfolder
			print "------------------------------------"
		
				
			#From verdict extract the paired TF-loc list								
			TFloc_list = parse_verdict(combined_verdict)
	
			#Prepare data for scatter plots
			data_array = link_data_in_array(TFloc_list, loc_dict, reg_dict)
			dataset_list.append(data_array)
			
	for i in range(0, len(dataset_list)):
		
		data_array_0 = dataset_list[0]
		data_array_1 = dataset_list[1]
		
		
	#Draw some scatter plots
	plot_title = "Predicted regulation from %s with cutoff %s"%(dataset, cutoff)
	scatterfilelocation = "%s/%s/scatter_combined_KS_%s_co%s.png"%(fa.mr_folder, fa.scatplot_folder, dataset, cutoff)
	print "------------------------------------"
	print "drawing %s"%scatterfilelocation
	print "------------------------------------"

	draw_combi_plot(data_array_0, data_array_1, scatterfilelocation, plot_title, chr_len_list, TF_chr_len_list)


if option == 5:
	

	
	for cutoff in cutoff_list:
			
		for dataset in exp_list:
			
			subfolder = "%s_%s/"%(dataset, cutoff)	
			combined_verdict = col.OrderedDict()
			print "------------------------------------"
			print "Analysing %s..."%subfolder
			print "------------------------------------"
			
			for chrom in chromosome:
				print "processing chromosome %s"%(chrom)
					
				filelocation = "%s/%s/%s/results2_%s_co%s_chr%s.txt"%(fa.main_folder, fa.enriched_folder, subfolder, dataset, cutoff, chrom)
				verdict = verify(flowerset, filelocation)
				#Combine dictionaries
				combined_verdict.update(verdict)
			
			print "------------------------------------"
			print "retrieving data for scatter plot..."
			print "------------------------------------"
		
				
			#From verdict extract the paired TF-loc list								
			TFloc_list = parse_verdict(combined_verdict)
	
			#Prepare data for scatter plots
			data_array = link_data_in_array(TFloc_list, loc_dict, reg_dict)
	
			#Draw some scatter plots
			plot_title = "Predicted regulation from %s with cutoff %s"%(dataset, cutoff)
			scatterfilelocation = "%s/%s/scatter_flower_%s_co%s.png"%(fa.mr_folder, fa.scatplot_folder, dataset, cutoff)
			print "------------------------------------"
			print "drawing %s"%scatterfilelocation
			print "------------------------------------"
			draw_plot(data_array, scatterfilelocation, plot_title, chr_len_list, TF_chr_len_list)
			
			
			
if option == 6:
	
	for cutoff in cutoff_list:
			
		for dataset in exp_list:
			
			subfolder = "%s_%s/"%(dataset, cutoff)	
			combined_verdict = col.OrderedDict()
			print "------------------------------------"
			print "Analysing %s..."%subfolder
			print "------------------------------------"
			
			for chrom in chromosome:
				print "processing chromosome %s"%(chrom)
					
				filelocation = "%s/%s/%s/results2_%s_co%s_chr%s.txt"%(fa.main_folder, fa.enriched_folder, subfolder, dataset, cutoff, chrom)
				verdict = verify(TF_fam_set, filelocation)
				#Combine dictionaries
				combined_verdict.update(verdict)
			
			print "------------------------------------"
			print "retrieving data for scatter plot..."
			print "------------------------------------"
				
			#From verdict extract the paired TF-loc list								
			TFloc_list = parse_verdict(combined_verdict)
	
			#Prepare data for scatter plots
			data_array1 = link_data_in_array(TFloc_list, loc_dict, reg_dict)
	
			
########################################################################
############################gxe#########################################
########################################################################
	for cutoff in cutoff_list:
			
		for dataset in exp_list:
			
			subfolder = "gxe_%s_%s/"%(dataset, cutoff)	
			combined_verdict = col.OrderedDict()
			print "------------------------------------"
			print "Analysing %s..."%subfolder
			print "------------------------------------"
			
			for chrom in chromosome:
				print "processing chromosome %s"%(chrom)
					
				filelocation = "%s/%s/%s/results_gxe_%s_co%s_chr%s.txt"%(fa.main_folder, fa.enriched_folder, subfolder, dataset, cutoff, chrom)
				verdict = verify(TF_fam_set, filelocation)
				#Combine dictionaries
				combined_verdict.update(verdict)
			
			print "------------------------------------"
			print "retrieving data for scatter plot..."
			print "------------------------------------"
				
			#From verdict extract the paired TF-loc list								
			TFloc_list = parse_verdict(combined_verdict)
	
			#Prepare data for scatter plots
			data_array2 = link_data_in_array(TFloc_list, loc_dict, reg_dict)
	
			#Draw some scatter plots
			plot_title = "normal vs gxe: Predicted regulation from %s with cutoff %s"%(dataset, cutoff)
			scatterfilelocation = "%s/%s/scatter_gxe_%s_co%s.png"%(fa.mr_folder, fa.scatplot_folder, dataset, cutoff)
			print "------------------------------------"
			print "drawing %s"%scatterfilelocation
			print "------------------------------------"
			draw_combi_plot(data_array1, data_array2, scatterfilelocation, plot_title, chr_len_list, TF_chr_len_list)














