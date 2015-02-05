import sys
import os
from sortedcontainers import SortedDict

from qtl.models import Gene

sys.path.append('/home/wouter/xenv_dj1.6/QTL')
#sys.path.append('/mnt/geninf15/prog/www/django/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

###############################################################################
#########################Create url Links######################################
###############################################################################




def GO_link_creator(GOterm, go_info_dict):
	"""
	This function will prepare hyperlinks for the identified GO terms
	to be displayed
	
	Also the name and definition of the GO term will be given
	
	"""
	
	main_url = "http://www.ebi.ac.uk/QuickGO/GTerm?id="
		
	link = main_url + GOterm
	
	if GOterm in go_info_dict:
		
		info = go_info_dict[GOterm]
		
	else:
		info = ['NA', 'NA']
	
	go_link_name_def_tuple = (link, GOterm, info[0], info[1])
		
	return go_link_name_def_tuple





def gene_link_creator(gene):
	"""
	Take a trait and create a link to a website
	
	Output is a list with link and gene name
	
	edit: changed to link to a google search because the TAIR website 
	began blocking the links.
	"""
	
	#main_url_begin = "http://www.arabidopsis.org/servlets/Search?type=general&search_action=detail&method=1&show_obsolete=F&name="
	#main_url_end = "&sub_type=gene&SEARCH_EXACT=4&SEARCH_CONTAINS=1"
	
	#url_start = "https://www.google.nl/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8#sourceid=chrome-psyapi2&ie=UTF-8&q="
	url_start = "https://www.google.nl/search?q="
		
	#link = main_url_begin + gene + main_url_end
	
	link = url_start + gene
	description = get_trait_description(gene)
	genelink_list = [link, gene, description]
	
	return genelink_list
	

	
###############################################################################	
#######################Prepare data for Display################################	
###############################################################################	


	
def get_trait_description(trait):
	"""
	Retrieve the functional description of the requested trait from the database
	"""
	
	try:
		trait_descr = Gene.objects.get(locus_identifier = trait).gene_model_description
		
	except Gene.DoesNotExist:
		trait_descr = "NA"
	
	return trait_descr	




def get_genes_with_go_from_qtl(gene_dict, go_gene_dict):
	"""
	Get the genes from each QTL for which a function was predicted
	
	In order to do this make a set() of the gene lists in go_dict
	And also make a set() of the gene list in qtl_gene_dict
	
	Compare the sets 
	Extend the gene lists
	The ones that are intersecting should remain
	
	This also has to be done for the gene list per go term
	
	gene_dict : dict[(trait, chr nmbr, marker start pos, marker end pos)] = [gene name, gene start, gene end] 
	go_gene_dict : dict[i, go, fu_p_value, adj(fu_p_val, go_frac_scA, go_frac_scB]) = [gene list]
	
	qtl_gogenes_dict : dict[(trait, chr nmbr, marker start pos, marker end pos)] = [gene list]
	"""
	
	qtl_gogenes_dict = SortedDict()
	
	for qtl in gene_dict:
		
		intersecting_genes = []
		gene_list111 = gene_dict[qtl]
		#some list comprehension, gotta love it
		#it takes all the gene names from [gene name, gene start, gene end]
		gene_list11 = [str(i[0]) for i in gene_list111]
		gene_list1 = set(gene_list11)
		
		for go in go_gene_dict:
			
			genelist2 = set(go_gene_dict[go])
			#the & signifies intersection
			#it is a special set() method
			temp_list = gene_list1 & genelist2
			
			intersecting_genes.extend(temp_list)
		
		#the intersecting genes are overlapping genes between genes from
		#QTLs with at least one enriched GO term and the entire genome
		qtl_gogenes_dict[qtl] = sorted(set(intersecting_genes))
		
	return qtl_gogenes_dict
		
	


def link_the_linked_links(gene_dict, go_gene_dict, go_definition_dict):
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
	
	#By merging gene_dict and go_gene_dict it is possible to arrange GO terms and enriched genes into eQTLS
	go_genes_qtl_dict : dict[(url link, GO term, name, definition)] = list of [link, gene, description]
	"""
	
	
	go_genes_qtl_dict = SortedDict()
	
	for go_info in go_gene_dict:
		
		go = go_info[1]
		go_fdr = go_info[3]
		
		#Create a list of url links for the GO terms
		#plus a GO definition
		#This will be the key in the dictionary
		link_go_tuple1 = GO_link_creator(go, go_definition_dict)
		link_go_tuple2 = (link_go_tuple1[0], link_go_tuple1[1], link_go_tuple1[2], link_go_tuple1[3],  go_fdr)

		intersecting_array = []
		go_gene_set = set(go_gene_dict[go_info])
		
		#print sorted(go_gene_set)
		for qtl in gene_dict:
			
			gene_list111 = gene_dict[qtl]
			#another list comprehension taking all the gene names..
			gene_list11 = [str(i[0]) for i in gene_list111]
			qtl_gene_set = set(gene_list11)
			#Another set() operation
			intersecting_genes = go_gene_set & qtl_gene_set
			
			intersecting_gene_link_list = []
			sorted_intersecting_genes = sorted(intersecting_genes)
			
			for gene in sorted_intersecting_genes:
				
				#create a list of url links for the gene names
				intersecting_gene_link_list.append(gene_link_creator(gene))
		
			intersecting_array.append(sorted(intersecting_gene_link_list))	
			
		go_genes_qtl_dict[link_go_tuple2] = intersecting_array
		

	#now we have a dictionairy containing for each go term a list of genes
	#each list of genes is separated into occurence in respective eQTLs
	
	return go_genes_qtl_dict
