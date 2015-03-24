import folder_assignments as fa

#Get TG-TF list from the enrichment lists
def get_TGTF_from_genelist(floc, TG_TF_ref, TG_list_ref, TF_list_ref):
	"""
	Extract data from several files
	The files are arranged per dataset and cutoff
	
	The data consists of traits with a genelist from their eQTL regions
	"""
	

	total_rel_major = []
	true_relations = []
	
	TF_set_ref = set(TF_list_ref)

													
	#These files contain the genelists from the eQTLs of every trait
	with open(floc, 'r') as fo:
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



def process_enrichment(trait_eqtl_genelist, TF_ref_set):
	"""
	"""
	TG_TF_relations = []
	for trait, eqtl, genelist in trait_eqtl_genelist:
		genelist_set = set(genelist)
		TF_in_genelist = TF_ref_set & genelist_set
		if TF_in_genelist:
			for TF in TF_in_genelist:
				TG_TF_relations.append([trait, TF])
				
	return TG_TF_relations



