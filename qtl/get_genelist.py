import sys
import os
import gc

from sortedcontainers import SortedDict

#The following lines enable this script to contact django
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD,Experiment,ExperimentMarker
from qtl.genelist_from_eqtl import (
	marker_logp_list, add_chromosome_and_position, retrieve_logp_above_cutoff, 
	check_region, iterate_marker_tuple, highest_tuple, set_region_from_adjacent_markers, get_chromosome_max, 
	check_regions_on_chromosome, combine_marker_region_data, find_genes_in_regions, normalize_object_data, 
	read_distinct_genes
	)
from qtl.go_enrichment import (
	read_once_csv, segregate_gene_list, annotate_from_csv, read_go_info_once_csv, 
	unique_GO_list, lets_see_tables, flatten_array, count_all_goterms, populate_contingency_table, 
	fish_for_python, extract_significant_result, correct_pvalues_for_multiple_testing, make_qtl_go_dict, 
	post_process_A, post_process_B)
	


#Main
main_folder = "/home/wouter/xenv_dj1.6/QTL"
script_folder = "validation_script"
result_folder = "validation"
mr_folder =  "%s/%s"%(main_folder, result_folder)
#Raw data
raw_folder = "raw_data"
filename_atreg = "%s/%s/AtRegNet.txt"%(result_folder, raw_folder)



	
def get_trait_name(trait):
	"""
	Retrieve the functional description of the requested trait from the database
	"""
	
	try:
		trait_descr = Gene.objects.get(locus_identifier = trait).primary_gene_symbol
		
	except Gene.DoesNotExist:
		trait_descr = "NA"
	
	return trait_descr	


def read_parse_atreg(fn):
	"""
	read the AtRegNet file into a dataobject
	"""
	with open(fn, 'r') as fo:
		data = fo.readlines()

	db = []
	for line in data:
		info = line.split('\t')
		
		nr = info[0].strip('"')
		gene_name = info[1].strip('"')
		TF = info[2].strip('"')
		family = info[3].strip('"')
		alt_name = info[4].strip('"')
		locus_id = info[5].strip('"')
		yesno = info[6].strip('"')
		direct = info[7].strip('"')
		confirmed = info[8].strip('"')
		biology = info[9].strip('"')
		activity = info[10].strip('"')
		article = info[11].strip('"')
		
		if len(TF.strip()) > 9:
			tf1 = TF[0:9]
			tf2 = TF[10:20]
			
			TF_split = [tf1, tf2]
			for tf_n in TF_split:
				cell = [nr,gene_name,tf_n,family,alt_name,locus_id,yesno,direct,confirmed,biology,activity,article]
				db.append(cell)
		
		else:		
		
			cell = [nr,gene_name,TF,family,alt_name,locus_id,yesno,direct,confirmed,biology,activity,article]
			db.append(cell)
		
	
	return db
		

def make_AtReg_dict(data):
	"""
	Make a dict with Tf as key and target as value
	"""
	
	AtReg_dict = SortedDict()
	for line in data:
		TF = line[2].upper()
		target = line[5].upper()
		
		if TF not in AtReg_dict:
			AtReg_dict[TF] = [target]
		else:
			temp_target = AtReg_dict[TF]
			temp_target.append(target)
			AtReg_dict[TF] = temp_target
		
	return AtReg_dict
		


def give_genelist_for_trait(trait, cutoff, exp, gxe_boolean):
	"""
	Check enriched genes for each trait for presence of regulator
	
	Perhaps write the genelists to a file for any threshold
	For later use...
	"""
	
	##########################################
	####Run functions towards enrichment######
	##########################################
	
	#Get genelist from eQTL
	l1 = marker_logp_list(trait, gxe_boolean, exp)
	l2 = add_chromosome_and_position(l1)
	l3 = retrieve_logp_above_cutoff(l2, cutoff)
	l4 = check_region(l3)
	l5 = iterate_marker_tuple(l4)
	l6 = set_region_from_adjacent_markers(l5, l2)	
	l7 = combine_marker_region_data(trait, l6)
	gene_dict = find_genes_in_regions(l7)

	#Count the number of found qtls
	qtls = len(gene_dict)
	
	#delete all objects after they become useless
	#because django is a memory horder
	#when gc.collect() is run it will delete all unreferenced dataobjects
	del l1
	del l2
	del l3
	del l4 
	del l5 
	del l6
	del l7
	
	print "trait %s has %s eQTLs"%(trait, qtls)
	
	if qtls != 0:
	
		gene_list = read_distinct_genes(gene_dict)

		return gene_list
		

	


def write_genelist_to_file(fn, gene, glist):
	"""
	
	"""
	print "writing file %s"%fn
	with open(fn, 'w') as fo:
		fo.write("trait: %s "%gene)
		fo.write("\n")
		loc_name = get_trait_name(gene)
		fo.write("trait name: %s "%loc_name)
		fo.write("\n")

		for g in glist:
			fo.write("TF: %s "%g)
			fo.write("\n")
			g_name = get_trait_name(g)
			fo.write("TF name: %s"%g_name)
			fo.write("\n")
	
	

def main():
	"""
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/get_genelist.py')		
	"""

		
	#Datasets
	exp_list = ['Ligterink_2014','Keurentjes_2007','Snoek_2012']

	#Variables
	#trait = "AT3G66652"
	dataset = exp_list[0]
	#chromosome = [1,2,3,4,5]
	cutoff = [3, 4.3, 6.7]
	coff = [cutoff[0]]
	gxe_bool = False	
	gxe = "False"
	
	loccer = []
	traitlist = ['AT1G62300', 'AT2G43010', 'AT3G27920', 'AT3G47640', 'AT4G36920']

	folder = "validation/gene_lists_AtRegNet/%s"%dataset
	
	if not os.path.exists(folder):
		os.mkdir(folder)

	TF_db = read_parse_atreg(filename_atreg)
	adict = make_AtReg_dict(TF_db)
	
	for co in coff:
		
		for loc in adict:
			gene_list = give_genelist_for_trait(loc, co, dataset, gxe_bool)
			TF_list = adict[loc]
			
			if gene_list:
				#Compare overlap between gene_list and TF_list
				gene_set = set(gene_list)
				TF_set = set(TF_list)
				
				overlap = gene_set & TF_set
				
				if overlap:
					loccer.append(loc)
					filename = "%s/genelist_%s_co%s_gxe%s_%s.txt"%(folder, dataset, co, gxe, loc)
					write_genelist_to_file(filename,loc, list(overlap))
					

	print loccer


main()


















					
