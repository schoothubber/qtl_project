import sys
import os
import gc

#The following lines enable this script to contact django
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

from qtl.models import Gene,Marker,LOD,Experiment,ExperimentMarker

from qtl.genelist_from_eqtl import marker_logp_list, add_chromosome_and_position, retrieve_logp_above_cutoff, check_region, iterate_marker_tuple, highest_tuple, set_region_from_adjacent_markers, get_chromosome_max, check_regions_on_chromosome, combine_marker_region_data, find_genes_in_regions, normalize_object_data, read_distinct_genes
from qtl.go_enrichment import read_once_csv, segregate_gene_list, annotate_from_csv, read_go_info_once_csv, unique_GO_list, lets_see_tables, flatten_array, count_all_goterms, populate_contingency_table, fish_for_python, extract_significant_result, multiple_test, make_qtl_go_dict, post_process_A, post_process_B
from qtl.prepare_for_display import get_genes_with_go_from_qtl

def trait_list_maker():
	"""
	locus_identifier = models.CharField(max_length=30,primary_key=True) #AT1G01480
    gene_model_name = models.CharField(max_length=30,blank = True)#AT1G01480.1
    gene_model_description = models.TextField(blank = True)
    gene_model_type = models.CharField(max_length = 40,blank = True)
    primary_gene_symbol = models.TextField(blank = True)
    all_gene_symbols = models.TextField(blank = True)
    chromosome = models.IntegerField(blank = True)
    start = models.IntegerField(blank = True)
    end = models.IntegerField(blank = True)
    orientation = models.BooleanField(blank = True) # Django True: sense strand False: anti-sense strand or MySQL 1: sense strand 0: antisense strand
	
	
	
	execfile('qtl/trait_list.py')
	"""

	exp_list = ['Ligterink_2014', 'Ligterink_2014_gxe', 'Keurentjes_2007', 'Snoek_2012']
	chromosome = [1,2,3,4,5]
	chrom = chromosome[1]
	
	for chrom in chromosome:
		traitlist = []
		for info in Gene.objects.all().filter(chromosome = chrom):
			trait = info.locus_identifier
			traitlist.append(trait)
			
		l = len(traitlist)
		s = len(set(traitlist))
		
		print "chr %s"%chrom
		print "l: %s"%l
		print "s: %s"%s

		
	#filename = "traitlist_chr%s.txt"%chrom
	#l = len(traitlist)
	#with open(filename, 'w') as fo:
		#for t in traitlist:
			#fo.write(t)
			#fo.write('\n')
		
	#print("traitlist for chromosome %s"%chrom)
	#print("number of traits: %s"%l)
		
		
trait_list_maker()
