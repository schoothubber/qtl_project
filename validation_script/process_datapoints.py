import folder_assignments as fa

#Get TG-TF list from the enrichment lists
def get_TGTF_from_genelist(dataset, cutoff, chromosome, TG_TF_ref, TG_list_ref, TF_list_ref):
	"""
	Extract data from several files
	The files are arranged per dataset and cutoff
	
	The data consists of traits with a genelist from their eQTL regions
	"""
	
	subfolder = "genelist_%s/"%dataset
	
	total_rel_major = []
	true_relations = []
	
	TF_set_ref = set(TF_list_ref)
			
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
					
					#Check if there are TFs in genelist
					gene_set = set(genelist)
					TF_in_genelist = gene_set & TF_set_ref
					
					if TF_in_genelist:
						for g_id_TF in TF_in_genelist:
							if [trait, g_id_TF] in TG_TF_ref:
								true_relations.append([trait, g_id_TF])
								
								
						total_rel_minor = count_total_relations(
														trait, TF_in_genelist, 
														TG_list_ref, TF_list_ref
														)
						total_rel_major.extend(total_rel_minor)


	return true_relations, total_rel_major


def count_total_relations(trait, TF_in_genelist, TG_list_ref, TF_list_ref):
	"""
	"""
	total_relations = []
	if TF_in_genelist:
		for g_id_TF in TF_in_genelist:				
				
			if trait in TG_list_ref and g_id_TF in TF_list_ref:
				total_relations.append([trait, g_id_TF])
				
	return total_relations



def process_enrichment(dataset, cutoff, chromosome, TF_set_ref, true_target_genes):
	"""
	This script should classify the datapoints into:
	-True Positives (purple)
	-False Positives (blue)
	"""


	enrichmentlist = []
	subfolder = "%s_%s/"%(dataset, cutoff)
	
	TG_TF_pred_major = []

	for chrom in chromosome:
		
		if dataset == 'Ligterink_2014_gxe':
			datset = 'Ligterink_2014'
			subfolder = "%s_%s/"%(datset, cutoff)

			filelocation = "%s/%s/v3_gxe_%s/results3_%s_co%s_chr%s.txt"%(
												fa.main_folder, fa.enriched_folder, 
												subfolder, datset, cutoff, chrom
												)
		
		else:
			filelocation = "%s/%s/v3_%s/results3_%s_co%s_chr%s.txt"%(
												fa.main_folder, fa.enriched_folder, 
												subfolder, dataset, cutoff, chrom
												)

		
		with open(filelocation, 'r') as fo:
			for line in fo:
				
				if line.startswith('trait'):
					
					if not enrichmentlist:
						
						trait = line[7:16]

					else:#if enrichmentlist

						TG_TF_pred_minor = identify_prediction(
													trait, enrichmentlist, 
													true_target_genes, 
													TF_set_ref
													)
						TG_TF_pred_major.extend(TG_TF_pred_minor)
						#reset the enrichmentlist for the next loop			
						enrichmentlist = []	
													
						#When the previous trait is processed
						#Copy the new trait
						trait = line[7:16]

					
				if line.startswith('AT'):
					enrichment_temp = line.split()
					#enrichmentlist will contain all enriched genes
					#Note that there are overlapping genes between 
					#various GO terms. Therefor turn the list into a
					#set to eliminate any redundancy
					enrichmentlist.extend(enrichment_temp)							

	
	return TG_TF_pred_major
	

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
