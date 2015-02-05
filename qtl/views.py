import os
import sys
import time
import csv
from collections import Counter
from sortedcontainers import SortedDict
from exceptions import EOFError

import scipy
from scipy import stats
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import StringIO

from django.shortcuts import render_to_response
from django.http import HttpResponse
from django.core.urlresolvers import reverse
from django.db import connection

from .models import Gene,Marker,LOD,Experiment,ExperimentMarker

from genelist_from_eqtl import marker_logp_list, add_chromosome_and_position, retrieve_logp_above_cutoff, check_region, iterate_marker_tuple, highest_tuple, set_region_from_adjacent_markers, get_chromosome_max, check_regions_on_chromosome, combine_marker_region_data, find_genes_in_regions, normalize_object_data, read_distinct_genes
from go_enrichment import read_once_csv, segregate_gene_list, annotate_from_csv, read_go_info_once_csv, unique_GO_list, lets_see_tables, flatten_array, count_all_goterms, populate_contingency_table, fish_for_python, extract_significant_result, correct_pvalues_for_multiple_testing, make_qtl_go_dict, post_process_A, post_process_B
from prepare_for_display import GO_link_creator, gene_link_creator, get_trait_description, get_genes_with_go_from_qtl, link_the_linked_links
from draw_DAG import OBOReader, GOTerm, GODag, after_colon, read_until




def BaseView(request):
	"""
	The view of the most basic homepage
	The URL http://127.0.0.1:8000/ will direct to wbase.html
	"""
	return render_to_response('wbase.html')



def SearchTraitView(request):
	"""
	When a trait and a cutoff value are given on the website, these values
	are retrieved from the database and fed to this view. 
	They are subsequently manipulated by other functions
	The output is sent to a nother template galled genelist.html
	
	The tic-toc was added to ascertain how long a search takes
	"""
	
	tic_t = time.clock() # starting point of total process
	
	#The template will show a checkbox, default will be unchecked
	#unchecked means 'gxe' is not in request.GET and thus we will go for
	#normale gene expression
	#checked means 'gxe' is in request.GET and thus we will go for 
	#environmental permutation
	if 'gxe' in request.GET:
		gxe_boolean = True
	if 'gxe' not in request.GET:
		gxe_boolean = False
		
	#The values 'trait_name' and 'cutoff_name' are given in the template
	if request.GET.get('trait_name') and request.GET.get('cutoff_name') and request.GET.get('experiment_name'):
		
		#Next check if there are any values given for the fisher test and
		#post processing procedures
		#If yes then use those values
		#If not then set them to the default values

		try:
			fisher_alpha_name_u = request.GET.get("fisher_alpha_name").strip()
			fisher_alpha_name = float(fisher_alpha_name_u)
			
		except:
			#default value:
			fisher_alpha_name = 0.05
		
		try:
			postA_name_u = request.GET.get("postA_name").strip()
			postA_name = float(postA_name_u)
			
		except:
			#default value:
			postA_name = 0.5
			
		try:
			postB_name_u = request.GET.get("postB_name").strip()
			postB_name = float(postB_name_u)
			
		except:
			#default value:
			postB_name = 0.01			
		
		
		tic = time.clock() # starting point of genelist retrieval
	
		trait_name_u = request.GET.get('trait_name').strip()
		trait_name = trait_name_u.encode("utf-8")
	
		cutoff_name_D = request.GET.get('cutoff_name').strip()
		cutoff_name = float(cutoff_name_D)
		
		experiment_name = request.GET.get('experiment_name').strip()
		#experiment_name = trait_name_u.encode("utf-8")

		trait_function = get_trait_description(trait_name)
		
		#Create url for GraphView
		#The variables for matplotlib will be supplied via the url
		Graph_url_base = request.build_absolute_uri(reverse('display_graph'))
		Graph_url_data = "?trait_name=%s&cutoff_name=%f&experiment_name=%s" %(trait_name, cutoff_name, experiment_name)
		Graph_url = Graph_url_base + Graph_url_data


		l1 = marker_logp_list(trait_name, gxe_boolean, experiment_name)
		l2 = add_chromosome_and_position(l1)
		l3 = retrieve_logp_above_cutoff(l2, cutoff_name)

		#If no markers are found above the chosen cutoff, than display this message
		if len(l3) == 0:
			message = "There are no markers with a LOD score higher than %f for trait %s" % (cutoff_name, trait_name)
			bottle = {
					"message" : message,
					"Graph_url" : Graph_url,
					}
			return render_to_response("wouter/display_individual_trait.html", bottle)
		
		else:
		
			l4 = check_region(l3)
			
			#this function uses the function 
			#highest_tuple
			l5 = iterate_marker_tuple(l4)
			
			#this function uses the functions 
			#get_chromosome_max and check_regions_on_chromosome
			l6 = set_region_from_adjacent_markers(l5, l2)
			
			l7 = combine_marker_region_data(trait_name, l6)
			
			#this function uses the function 
			#normalize_object_data
			gene_dict = find_genes_in_regions(l7)
			
			gene_list = read_distinct_genes(gene_dict)
			

			
			#Place the output into a dictionary
			#Give each key the same name as the variable for convenience
			#The dictionary will be passed to the template
			#Inside the template the key's are used to display the
			#corresponding values
			
			toc = time.clock()
			print "Time to retrieve a gene list is %f seconds" % (toc-tic)
			
			
			
			##############################################################
			#####################annotate genes###########################
			##############################################################
	
			
			#get the location of the GO annotation file
			filename_annotate = 'At_geneGOlists.csv'
			module_dir = os.path.dirname(__file__)  # get current directory
			file_path_annotate = os.path.join(module_dir, 'documents', filename_annotate)
			
			#read the rows from the csv file mentioned above in filename
			data_dict = read_once_csv(file_path_annotate)
			
			absent_genes, present_genes = segregate_gene_list(data_dict, gene_list)


			print "number of genes in gene_list = %d" %len(gene_list)
			print "number of genes in data_dict = %d" %(len(data_dict) - len(gene_list))	
			
			gene_inside_qtl_dict, gene_outside_qtl_dict, tot_in_qtl, tot_out_qtl = annotate_from_csv(data_dict, present_genes)
			print "missings = %d" % len(absent_genes)
			
			print "number of genes in qtl = %d and %d" %(len(gene_inside_qtl_dict), len(present_genes))
			print "number of genes out qtl = %d" %len(gene_outside_qtl_dict)
		
			golist_unique = unique_GO_list(gene_inside_qtl_dict)
			
			
			#get the location of the GO annotation file
			filename_bio = 'go_bio_process.csv'
			file_path_bio = os.path.join(module_dir, 'documents', filename_bio)
			#This reads all Biological Process GO terms with descriptions	
			go_info_dict = read_go_info_once_csv(file_path_bio)
			
			
			

			##############################################################
			##############################################################		


			#Count all genes with and without GO annotations
			#Count all genes inside and outside the QTLs
			
			#Prepare for counted genes inside qtl:
			golist_in_flat = flatten_array(gene_inside_qtl_dict)
			c_go_in_qtl_dict = count_all_goterms(golist_in_flat)
			
			#Prepare for counted genes outside qtl:
			golist_out_flat = flatten_array(gene_outside_qtl_dict)
			c_go_out_qtl_dict = count_all_goterms(golist_out_flat)
			
			#Create contingency tables for the Fishers exact test						
			tic = time.clock()
			c_array, total_genes = populate_contingency_table(c_go_in_qtl_dict, c_go_out_qtl_dict, golist_unique, tot_in_qtl, tot_out_qtl)
			toc = time.clock()
			print "Time to make contingency table is %f seconds" % (toc-tic)			
			
				
				
			##############################################################
			####################Fishing with SciPy########################
			##############################################################
			
			
			
			#perform the fisher exact test on all created contingency
			#tables
			fisher_python = fish_for_python(c_array)
			
			#only allow the results with p values below fisher_alpha to pass
			#fisher_alpha default is 0.05 unless another value is given
			#front end
			significant_info = extract_significant_result(fisher_python, fisher_alpha_name)
			
			#perform a multiple test
			enriched_golist = correct_pvalues_for_multiple_testing(significant_info)
			
			#Post processing
			qtl_go_dict = make_qtl_go_dict(gene_inside_qtl_dict, gene_dict)
					

			#A
			tic = time.clock()
			approved_golistA = post_process_A(enriched_golist, qtl_go_dict, postA_name)
			toc = time.clock()
			print "Time to approve enrichment A is %f seconds" % (toc-tic)
			
			
			#gene_list_for_specific_go = get_gene_with_go_from_qtl(approved_golistA, gene_inside_qtl_dict, gene_outside_qtl_dict)
			
			
			#B
			tic = time.clock()
			godict_full = post_process_B(approved_golistA, gene_inside_qtl_dict, gene_outside_qtl_dict, postB_name)
			toc = time.clock()
			print "Time to approve enrichment B is %f seconds" % (toc-tic)
			print "length of godict = %d" %len(godict_full)
			
			#implement function that compares go terms and arranges according to parent child
			#relationships if possible
			
			go_term = "GO:0061024"
			#Create url for DAGView
			#The variables for graphviz will be supplied via the url
			DAG_url_base = request.build_absolute_uri(reverse('display_godag'))
			DAG_url_data = "?goterm_name=%s" %(go_term)
			DAG_url = DAG_url_base + DAG_url_data
		
		
			#GO links:
			link_for_trait = gene_link_creator(trait_name)
			
			
			#go_gene_link_dict = link_the_linked_links2(godict_reduced)
			
			qtl_gogenes_dict = get_genes_with_go_from_qtl(gene_dict, godict_full)
			
			go_genes_qtl_dict = link_the_linked_links(gene_dict, godict_full, go_info_dict)
			
			
			
			#c_array_display = lets_see_tables(c_array)
			
			

					
			stats_dict = {
						"link_for_trait" : link_for_trait,
						"Graph_url" : Graph_url,
						"DAG_url" : DAG_url, 
						"cutoff_name" : cutoff_name,
						"qtl_gogenes_dict" : qtl_gogenes_dict,
						"go_genes_qtl_dict" : go_genes_qtl_dict,
						}
			
			
			toc_t = time.clock()
			print "Time to perform entire process is %f seconds" % (toc_t-tic_t)
			


			
			#return display_meta(request)
			return render_to_response("wouter/display_individual_trait.html", stats_dict)

		
	else:
		exps = Experiment.objects.all().values_list('experiment_name',flat=True)
		return render_to_response('wouter/search_individual_trait.html',{'exps':exps})
		#return render_to_response('wouter/search_individual_trait.html')







def MultipleTraitView(request):
	"""
	This view will display the amount of found qtl and GO enrichments
	for each trait for a given cutoff value
	
	First: check number of qlt, should be bigger than 1
	
	Second: check number of found GO enrichments
	
	"""
	
	#return display_meta(request)
	
	if 'gxe' in request.GET:
		gxe_boolean = True
	if 'gxe' not in request.GET:
		gxe_boolean = False
		
		
	if request.GET.get('trait_name') and request.GET.get('cutoff_name'):
		tic = time.clock()
		
		trait_query = request.GET.get('trait_name')
		trait_list1 = trait_query.split()
		trait_list2 = []
		
		for trait in trait_list1:
			trait.strip(",")
			trait.strip("'")
			trait_list2.append(trait)
		
		combined_stats = []
		
		#translate cuttoff
		cutoff_name_D = request.GET.get('cutoff_name')
		cutoff_name = float(cutoff_name_D)
		
		experiment_name = request.GET.get('experiment_name').strip()
		
		
		#get list of all traits
		#traits = Gene.objects.values_list('locus_identifier', flat=True).order_by('locus_identifier')

		#find number of qtl for trait
		for trait in trait_list2:
			
			trait_name = trait.encode("utf-8")
			
			#First search for the amount of QTL's
			l1 = marker_logp_list(trait_name, gxe_boolean, experiment_name)
			l2 = add_chromosome_and_position(l1)
			l3 = retrieve_logp_above_cutoff(l2, cutoff_name)
			l4 = check_region(l3)
			l5 = iterate_marker_tuple(l4)
			l6 = set_region_from_adjacent_markers(l5, l2)
			l7 = combine_marker_region_data(trait_name, l6)
			
			#Howmany qtl for this trait?....
			number_of_qtl = len(l7)
			
			#Second search for the amount of GO enrichments
			gene_dict = find_genes_in_regions(l7)
			gene_list = read_distinct_genes(gene_dict)
			
			number_of_genes = len(gene_list)
			
			
			combined_stats.append([trait_name, number_of_qtl, number_of_genes])
			
			
			#print "trait %s has %d qtls with %d genes" % (trait_name, number_of_qtl, number_of_genes)

		toc = time.clock()
		print "Time to perform trait check is %f seconds" % (toc-tic)
		
		
		
		stats_dict = {
					"cutoff_name" : cutoff_name,
					"combined_stats" : combined_stats
					}
	
		
		return render_to_response("wouter/display_multiple_traits.html", stats_dict)
	
	
	else:
		exps = Experiment.objects.all().values_list('experiment_name',flat=True)
		return render_to_response('wouter/search_multiple_trait.html',{'exps':exps})
		#return render_to_response("wouter/search_multiple_trait.html")





def Graphview(request):
	"""
	A graph of the LOD scores versus the markers
	
	The graph is activated through the SearchTraitView via its namespace
	
	The Graphview gets its variables via the url adress
	
	"""
	tic = time.clock()
	if request.GET.get('trait_name') and request.GET.get('cutoff_name') and request.GET.get('experiment_name'):
		
		trait_name_u = request.GET.get('trait_name')
		trait_name = trait_name_u.encode("utf-8")
	
		cutoff_name_D = request.GET.get('cutoff_name')
		cutoff_name = float(cutoff_name_D)
		gxe_boolean = False
		
		experiment_name = request.GET.get('experiment_name').strip()
		
		data = marker_logp_list(trait_name, gxe_boolean, experiment_name)
		
		#labels = [str(i[1]) for i in data]
		x = [int(i[0]) for i in data]
		y = [float(i[2]) for i in data]
		
		fig = Figure(figsize = (15,5))
		axis = fig.add_subplot(1,1,1)
		axis.set_title("LOD scores for trait %s with cutoff %s" %(trait_name, cutoff_name))
		axis.set_xlabel("Marker")
		axis.set_ylabel("LOD")
		axis.grid(True)
		axis.axhline(0, color = 'k')
		axis.axhline(cutoff_name, color = 'r')
		axis.axhline(-cutoff_name, color = 'r')
		axis.autoscale(enable = True)
		axis.plot(x, y)
		
		canvas = FigureCanvas(fig)
		
		imgData = StringIO.StringIO()
		canvas.print_png(imgData)
		response = HttpResponse(imgData.getvalue(), content_type = 'image/png')
		
		toc = time.clock()

		print "time for graph construction is %f" %(toc-tic)
		return response
			

	else:
		exps = Experiment.objects.all().values_list('experiment_name',flat=True)
		return render_to_response('wouter/search_graph_trait.html',{'exps':exps})
		#return render_to_response("wouter/search_graph_trait.html")
	


def DAGView(request):
	"""
	Construct a Directed Acyclic Graph around one GO term
	"""
	
	if request.GET.get('goterm_name'):
		goterm_name_u = request.GET.get('goterm_name')
		goterm_name = goterm_name_u.encode("utf-8")
		
		filename_obo = 'go-basic.obo'
		module_dir = os.path.dirname(__file__)  # get current directory
		file_path_obo = os.path.join(module_dir, 'documents', filename_obo)
				
		g = GODag(file_path_obo)
		
		rec = g.query_term(goterm_name, verbose=True)
		#imgData = StringIO.StringIO()
		s = g.draw_lineage([rec], gml=False, draw_parents=True, draw_children=True)
		
		response = HttpResponse(s , mime_type = 'image/svg')

		return response


def StoreCsvView(request):
	"""
	Work in progress...
	"""
	
	# Create the HttpResponse object with the appropriate CSV header.
	response = HttpResponse(content_type='text/csv')
	response['Content-Disposition'] = 'attachment; filename="somefilename.csv"'
	
	writer = csv.writer(response)
	
	#Iterate through all lists/dictionaries/tuples, anything containing
	#information to write to the csv file
	writer.writerow(['First row', 'Foo', 'Bar', 'Baz'])
	writer.writerow(['Second row', 'A', 'B', 'C', '"Testing"', "Here's a quote"])
	
	return response

# AT1G03530
###############################################################################



def display_meta(request):
	"""
	Display all the data (keys and values) that is inside the request
	just for kicks
	"""
	values = request.META.items()
	values.sort()
	html = []
	for k, v in values:
		html.append('<tr><td>%s</td><td>%s</td></tr>' % (k, v))
	return HttpResponse('<table>%s</table>' % '\n'.join(html))
		
	
		
###############################################################################
############################TESTING VIEW#######################################
###############################################################################
###############################################################################
		

def OutputDataView(request):
	"""
	When a trait and a cutoff value are given on the website, these values
	are retrieved from the database and fed to this view. 
	They are subsequently manipulated by other functions
	The output is sent to a nother template galled genelist.html
	
	The tic-toc was added to ascertain how long a search takes
	"""
	
	#The template will show a checkbox, default will be unchecked
	#unchecked means 'gxe' is not in request.GET and thus we will go for
	#normale gene expression
	#checked means 'gxe' is in request.GET and thus we will go for 
	#environmental permutation

		
	if 'gxe' in request.GET:
		gxe_boolean = True
	if 'gxe' not in request.GET:
		gxe_boolean = False
		
	#The values 'trait_name' and 'cutoff_name' are given in the template
	if request.GET.get('trait_name') and request.GET.get('cutoff_name') and request.GET.get('experiment_name'):
		
	
		trait_name_u = request.GET.get('trait_name').strip()
		trait_name = trait_name_u.encode("utf-8")
	
		cutoff_name_D = request.GET.get('cutoff_name').strip()
		cutoff_name = float(cutoff_name_D)
		
		experiment_name = request.GET.get('experiment_name').strip()
	
		trait_function = get_trait_description(trait_name)

		l1 = marker_logp_list(trait_name, gxe_boolean, experiment_name)
		l2 = add_chromosome_and_position(l1)
		l3 = retrieve_logp_above_cutoff(l2, cutoff_name)

		#If no markers are found above the chosen cutoff, than display this message
		if len(l3) == 0:
			message = "There are no markers with a LOD score higher than %f for trait %s" % (cutoff_name, trait_name)
			bottle = {"message" : message}
			return render_to_response("wouter/display_individual_trait.html", bottle)
		
		else:
		
			l4 = check_region(l3)
			
			#this function uses the function 
			#highest_tuple
			l5 = iterate_marker_tuple(l4)
			
			#this function uses the functions 
			#get_chromosome_max and check_regions_on_chromosome
			l6 = set_region_from_adjacent_markers(l5, l2)
			
			#l7 = combine_marker_region_data(trait_name, l6)
			
			#this function uses the function 
			#normalize_object_data
			#gene_dict = find_genes_in_regions(l7)
			
			#gene_list = read_distinct_genes(gene_dict)
			
			
			
			######################
			markers = Marker.objects.values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
			cm_data = []
			for marker in markers:
		
				cm_m = Marker.objects.get(marker_name = marker).marker_cm
				cm_data.append(cm_m)
				
			######################
			
			
			
			gene_dict = {
						"trait_name" : trait_name,
						"cutoff_name" : cutoff_name,
						"l2" : l2,
						"l3" : l3,
						"l4" : l4,
						"l5" : l5,
						"l6" : l6,
						#"gene_dict" : gene_dict,
						#"gene_list" : gene_list,
						"cm_data" : cm_data,
						}



			return render_to_response("wouter/display_output_data.html", gene_dict)

		
	else:
		exps = Experiment.objects.all().values_list('experiment_name',flat=True)
		return render_to_response('wouter/search_graph_trait.html',{'exps':exps})
		#return render_to_response('wouter/search_individual_trait.html')
		
		
		
	
