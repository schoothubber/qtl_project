import os
import time
import csv
from collections import Counter
from sortedcontainers import SortedList, SortedSet, SortedDict
import random

import rpy2.robjects as ro
import numpy as np
import scipy
from scipy import stats
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import PIL #Python Imaging Library
from PIL import Image
import StringIO

from django.shortcuts import render_to_response
from django.http import HttpResponse
from django.core.urlresolvers import reverse
from django.db import connection

from .models import Gene,Marker,LOD,Experiment,ExperimentMarker



###############################################################################
#################################WOUTERS#######################################
###################################CODE########################################
###############################################################################

class Query:
	
	
	def __init__(self, trait, gxebool, experiment, cutoff, fisher_alpha, qtl_fraction, function_fraction):
		"""
		Initialize a query
		"""
		self.trait = trait
		self.gxebool = gxebool
		self.experiment = experiment
		self.cutoff = cutoff
		self.fisher_alpha = fisher_alpha
		self.qtl_fraction = qtl_fraction
		self.function_fraction = function_fraction
	
	
	def marker_logp_list(self):
		'''
		Create a list of tuples containing markers with corresponding -log p value.
		Both are subscribed to the requested trait.
		The markers and -log p values are ordered based on the markers' location on the
		chromosomes and their genetic mapping in centiMorgan.
		 
		The list is enumerated, meaning that an integer is added to
		the tuple. 
		
		The enumeration will eventually help to identify marker regions
		in the check_if_markers_consecutive function
		'''
		
		#First get all the marker names from the database and order them
		#according chromosomal position and genetic mapping

		markers_full = ExperimentMarker.objects.filter(experiment_name = self.experiment)
		markers = markers_full.values_list('marker_name', flat = True)
		
		#print connection.queries
		
		marker_logp_tuple = ()
		tuple_list = []
		 
		# iterate of the list of markers and get the -log p value for each marker
		for i in range(0, len(markers)):
			#print markers[i]
			
			try:
				#Because there is only one object that will match the query
				#the get() method is used on a Manager (in this case objects)
				#to return the object directly
				#gxe differentiates between normal gene expression (False)
				#and environmental permutation (True)
				logp1 = LOD.objects.get(locus_identifier = self.trait, marker_name = markers[i], gxe = self.gxebool, experiment_name = self.experiment)

				#From the returned object extract the value from the attribute LOD_score
				logp2 = logp1.LOD_score
				
				#Change the type from "Class decimal.Decimal" to float
				logp3 = float(logp2)#the actual -log p value
				
				#add the i for enumeration
				#perform encode("utf-8") on each marker to reduce the size
				marker_logp_tuple = (i, markers[i].encode("utf-8"), logp3)
				tuple_list.append(marker_logp_tuple)
				
			except LOD.DoesNotExist:
				pass
		
		return tuple_list
			#trait, [(int, marker name, -logp)]

		
	def add_chromosome_and_position(self, tuple_list): #[(int, marker, -logp)]
		"""
		To correlate the marker positions on the chromosome we need the 
		markers' chromosome number and physical position on that chromosome
		These extra values are taken from the Marker table, using the marker_name
		as Foreign Key
		"""
		#these values are added here already so that later functions can be
		#less complicated
		
		info_tuple = ()
		info_list = []
		
		#Iterate through the list with tuples
		for data_tuple in tuple_list:
			
			chromosome_number = None
			physical_pos = None
					
			#for each tuple add the chromosome number and the physical position
			#these are taken from the Marker table with marker_name as foreign key
			chromosome_number = Marker.objects.get(marker_name = data_tuple[1]).marker_chromosome
			physical_pos = Marker.objects.get(marker_name = data_tuple[1]).marker_phys_pos
			
			#The physical position is given as a type class decimal.Decimal with digit(s) in front of 
			#the comma and six of the last digits behind the comma plus some zero's
			#multiplying by a million will give the physical location in basepairs
			#Conversion to an integer is necessary to change the type
			pp = int((physical_pos) * 1000000)
			
			#so that each tuple will become:
			#(integer, marker, -logp, chromosome number, physical position)
			info_tuple = (data_tuple[0], data_tuple[1], data_tuple[2], int(chromosome_number), (pp))
			info_list.append(info_tuple)
		
		return info_list #list of tuples		


	def retrieve_logp_above_cutoff(self, tuple_list):
		'''
		Take all the Markers and -logp values above a certain value,
		(Usually a cutoff of 2.3 is seen as significant)
		
		The cutoff is given as input in the app
		'''

		cutoff_list = []
		# Iterate through the tuple_list and exclude all tuples with a
		#-logp value lower than the cutoff
		for info_tuple in tuple_list:#for each tuple (int, marker, lodscore)
			
			
			#Use abs() because there are negative -logp values present
			#The minus signs were added to discriminate between the parent alleles
			# + indicates Bay and - indicates Sha
			if abs(info_tuple[2]) >= self.cutoff:		
					
				cutoff_list.append(info_tuple)
			else:
				continue
				
		return cutoff_list #(int, marker, -logp)
	

	def check_region(self, cutoff_list):
		#(integer, marker, -logp, chromosome number, physical position)
		"""
		Take a list of tuples.
		The integers in the list are checked whether they increase in 
		consecutive order (with an increment of 1).
		Tuples that group together in this way indicate a region of markers
		on any chromosome.
		"""


		#LOGICAL STEPS TO BE TAKEN:
		#check each value with its consecutive value
		#if difference is 1 put both values in a group_list
		#if the difference is no longer 1, it is no longer consecutive
		#i.e. no more linkage between markers
		#put group_list in a group_array
		#start new group_list with the new value
		#this also happens when the chromosome number changes between two
		#consecutive markers
		#and continue
		
		group_list = []
		group_array = []
		
		#Iterate over all the tuples in the cutoff_list
		for i in range(0, len(cutoff_list)):
			
			#the last value in the list cannot be used here
			#hence the -1
			#or else an Indexerror will occur!
			if i != len(cutoff_list)-1:
				
				#check each value with its consecutive value
				#if difference is 1 put integer[i] in a group_list
				#AND check if the chromosome numbers are equal
				#because the region has to be on one singular chromosome
				if (cutoff_list[i][0] == cutoff_list[i+1][0]-1) and ((cutoff_list[i][3]) == (cutoff_list[i+1][3])):			
					group_list.append(cutoff_list[i])
				   
				else:#if the difference is no longer+1 it is no longer a region
					#store the tuple together in a list
					group_list.append(cutoff_list[i])
					#put group_list in a group_array
					group_array.append(group_list)
					#And reset the group_list
					group_list = []

			else:#to prevent an Indexerror the last value is handled here
				#check if last one is 1 more than the second last one
				if (cutoff_list[i][0] == cutoff_list[i-1][0]+1) and ((cutoff_list[i][3]) == (cutoff_list[i-1][3])):
					group_list.append(cutoff_list[i])
					#add last tuple to the group_list
					group_array.append(group_list)
					#finished! dump the last group_list in the group_array

				else:#if the difference is not 1
					group_list = [cutoff_list[i]]
					#start a new group_list containing just the last integer
					group_array.append(group_list)#add it to the other groups

		return group_array
   

	def iterate_marker_tuple(self, marker_region_list):
		"""
		Iterate over a list containing lists with either one or more
		tuples. 
		
		Also a single tuple is sent, however, it will be identified quickly
		as the maximum in the group.(Considering its the only one!)
		
		A list of 2 or more tuples will be fed to the function highest_tuple. 
		"""
		
		maximum_markers = []
		
		#Iterate over the lists containing the tuple(s)
		for region in marker_region_list:
		
			maximum_markers.append(self.highest_tuple(region))
		
		return maximum_markers

	
	def highest_tuple(self, region):
		"""
		region = list of tuples
		tuple = (integer, marker name, -logp value, marker position)
		"""
		#Immediately set the first (or only) tuple in the list as maximum
		max_tuple = region[0]
		
		#If the region consist of one marker tuple then the work of this function
		#is easy, return the tuple immediately!
		if len(region) == 1:
			
			return max_tuple
			
		else:
			#If there are multiple tuples:
			#Iterate through the marker tuples in the list and compare 
			#the size of the -log p value of each marker tuple
			for i in range(0, len(region)-1):
				
				#use the abs() function because the -logp values can be negative
				if abs(region[i][2]) < abs(region[i+1][2]):
					
					max_tuple = region[i+1]
					
				else:
					
					max_tuple = region[i]
					
			return max_tuple
    
        
	def set_region_from_adjacent_markers(self, maximum_list, marker_list):
		"""
		The adjacent markers will be the actual region.
		
		Exceptions are made for cases where the highest scoring marker is
		the first one on the first chromosome or the last marker on the last
		chromosome.
		
		In case an adjacent marker lies on a different chromosome,
		the region border should be either the beginning or the end of
		the markers' chromosome
		
		Note that if the highest scoring marker lies at the beginning or the end
		of a chromosome, the neighbouring location cannot be another marker,
		but instead it should be the beginning or end of the chromosome.
		This function does not take this into account.
		But it is rectified in the function check_regions_on_chromosome via 
		this function
		"""
		
		#The eventual region begins and ends will be appended to this list
		actual_region_list = []
		
		#Iterate over the tuples in the list
		#each tuple: (integer, marker, -logp, chromosome, physical location)
		#Maintain datastructure for easy use in later functions
		for i in range(0, len(maximum_list)):
			
			region_start = None
			region_end = None
			
			#Check if the marker is located at the start of chromosome 1
			if maximum_list[i][0] == 0:
				#region = marker[+1] and position 0 on chr 1
				region_start = (None, None, None, maximum_list[i][3], 0)
				
				#This variable becomes an integer that is one higher than the
				#integer of the marker tuple that had the highest -log p value
				index_end = maximum_list[i][0]+1
				region_end = marker_list[index_end]
				
				actual_region_list.append([region_start, region_end])
			
			#Check if the marker is located at the end of the last chromosome
			elif maximum_list[i][0] == len(marker_list)-1:
				#region is the marker tuple with index_number minus one
				#untill the position chr_max
				chr_max = self.get_chromosome_max(maximum_list[i])
				
				#This variable becomes an integer that is one less than the
				#integer of the marker tuple that had the highest -log p value
				index_start = maximum_list[i][0]-1
				
				#take back the marker tuple from the original list
				region_start = marker_list[index_start]
				
				#set the new tuple for the end of a chromosome
				#where chr_max is the highest physical position for a gene
				#on a certain chromosome
				region_end = (None, None, None, maximum_list[i][3], chr_max)
				
				
				actual_region_list.append([region_start, region_end])
				
			else:#If we arrive here then the marker borders are not at the
				#edges of a chromosome
				corr_start, corr_end = None, None
				
				index_start = maximum_list[i][0]-1
				region_start = marker_list[index_start]
				
				
				index_end = maximum_list[i][0]+1
				region_end = marker_list[index_end]
				
				#Give the region borders to check_regions_on_chromosome.
				#The output is either the same data or an adjusted version.
				#This depends on the chromosomal positions of the borders.
				#However the datastructure of the tuple is kept
				corr_start, corr_end = self.check_regions_on_chromosome(region_start, region_end, marker_list)
				
				actual_region_list.append([corr_start, corr_end])
		
				
		return actual_region_list


	def get_chromosome_max(self, marker_tuple):
		"""
		The tuple being fed is always from the last marker on any
		given chromosome.
		
		Get the highest physical position on the chromosome to act as
		the highest region border
		"""
		#input(integer, marker, -logp, chromosome, physical location)
		#marker_tuple[3] is the chromosome number
		#All the end positions of the genes from a certain chromosome are 
		#selected, subsequently the highest value from that list is taken
		chr_size_list = Gene.objects.values_list("end", flat=True)
		chr_spec_list = chr_size_list.filter(chromosome = marker_tuple[3])
		chr_max = max(chr_spec_list)
		chr_max_int = int(chr_max)
		
		return chr_max_int
		

	def check_regions_on_chromosome(self, region_start, region_end, marker_list):
		#(integer, marker, -logp, chromosome, physical location)
		"""
		Regions cannot contain markers from different chromosomes
		This function will compare the chromosomal positions of two
		markers that were chosen as a region.
		
		If the chromosome numbers between the two tuples do not match,
		make the appropriate corrections:
		Depending on the location of the highest scoring marker and the
		marker on the next or previous chromosome, the start or end location
		of a region will be set to the start or end of the chromosome of the 
		highest scoring marker.
		"""
		
		
		if region_start[3] == region_end[3]:
			#perfect! The borders are on the same chromosome
			#just return the data and do nothing else
			return region_start, region_end
			
		#If the borders are on different chromosomes however,
		#revert to the first tuple containing the highest -logp value
		else:
			revert_index = region_start[0]+1
			#Note that is doesnt matter here if region_start[0]+1 is used
			#or region_end[0]-1. The answer would be the same integer
			revert_tuple = marker_list[revert_index]
			
			#if the topmarkers chromosome number is equal to the start marker
			#then keep the start marker as it is and alter the end marker
			if revert_tuple[3] == region_start[3]:
				corr_start = region_start
				#The region end should be on the same chromosome as the
				#region start, thats why region_start is fed to this function
				chr_max = self.get_chromosome_max(region_start)
				corr_end = (None, None, None, region_start[3], chr_max)
				return corr_start, corr_end
			
			#if the topmarkers chromosome number is equal to the end marker
			#then keep the end marker as it is and alter the start marker
			elif revert_tuple[3] == region_end[3]:
				#The region start should be on the same chromosome as the
				#region end.
				#The None values should not matter, this new tuple only has
				#a physical position on a chromosome and is not related to
				#any marker and thus does not have a LOD score
				corr_start = (None, None, None, region_end[3], 0)
				corr_end = region_end
				return corr_start, corr_end
			

	def combine_marker_region_data(self, region_tuple_list):
		"""
		In order to resume with the Go Enrichment it is easier to combine
		each region (consisting of two tuples at this moment) to one tuple
		with the following structure:
		(trait, chromosome number, marker start position, marker end position)
		This final tuple will be used as a dictionary key.
		The corresponding dictionary value will be a list of genes 
		from that particular marker region
		"""
		combined_data = []
		#region(integer, marker, -logp, chromosome, physical location)
		for region in region_tuple_list:
			combined_tuple = (self.trait, region[0][3], region[0][4], region[1][4])
			combined_data.append(combined_tuple)
			
			
		#combined_data(trait, chromosome , start physical location, end physical location)
		return combined_data
		
		
	def find_genes_in_regions(self, region_list):
		"""
		The region_list will contain tuples with the following information:
		(trait, chromosome number, marker start position, marker end position)
		
		This information is used to collect lists of genes that are positioned
		on each region.(between the 2 physical start and end positions)
		[gene name, gene start, gene end]
		
		The tuple and gene lists will be placed together in a dictionary with
		the tuple as key and the list as value
		"""
		
		gene_dict = SortedDict()
		
		complete_gene_list = Gene.objects.values_list("locus_identifier","chromosome", "start", "end")
		
		for region in region_list:
			#Only look at the genes on a specific chromosome
			#The chromosome number is indicated on region[1]
			chromosome_spec_list = complete_gene_list.filter(chromosome = region[1])
			#Exclude all genes with a physical location in front of the start 
			#position
			remove_genes_ahead = chromosome_spec_list.exclude(start__lt = region[2])
			#Exclude all genes with a physical location behind the end
			#position
			remove_genes_behind = remove_genes_ahead.exclude(end__gt = region[3])
			
			#normalize the database data, meaning that we use the data as
			#basic python types instead of the django class types
			normal_gene_list = self.normalize_object_data(remove_genes_behind)
			
			#add each region and list of genes to the dictionary gene_dict
			#with the tuple as a key and the gene list as a value
			gene_dict[region] = normal_gene_list
			
		return gene_dict
	

	def normalize_object_data(self, data_list):
		"""
		The extracted data from the database needs to be parsed 
		data_list has the structure 
		[gene name,chromosome number, gene start, gene end]
		With data types that were derived from the model classes
		The output should become 
		[gene name, gene start, gene end]
		With normalized data types
		"""	
		normal_gene_list = []
		
		for data in data_list:
			normal_gene_list.append([data[0].encode("utf-8"), int(data[2]), int(data[3])])
		
		return normal_gene_list


###############################################################################
###############################################################################


def read_distinct_genes(genedict):
	"""
	All genes from all regions will be placed in one big list.
	"""
	
	#read all the genes from the found regions into a new list
	#This list contains only genes which have to be unique
	#i.e. no gene can be mentioned more than once
	gene_list = []
	
	for qtl in genedict:
		
		for genedata in genedict[qtl]:
			
			if genedata[0] not in gene_list:
				gene_list.append(genedata[0])
			
			
	return gene_list




def read_once_csv(filename):
	"""
	Read a csv file and place the contents in a dictionary
	
	The csv file contains rows of data
	In each row the first index holds a gene name
	The rest of the row holds the GO annotations
	"""
	
	
	data_dict = {}
	
	with open(filename, 'rb') as csvobject:
		reader = csv.reader(csvobject)
		for row in reader:
			
			#key is the gene name
			key = row[0]
			#value is a list of go terms
			value = row[1:]
			
			data_dict[key] = value
	
	
	return data_dict



def segregate_gene_list(data_dict, gene_list):
	"""
	Take a dictionary, containing all the genes for which 1 or more GO
	terms are known, derived from a csv file.
	
	And take a list of genes that were located within the identified QTL 
	regions.
	
	Using pythons set() mechanics segregate the list of genes into 2 new
	lists of genes:
	1: a list of genes present in the data file
	2: a list of genes absent in the data file
	
	If some genes are not mentioned in the csv file, there is no point in
	using them anymore. This will make further processes shorter/faster
	"""

	
	genes_from_file = set(data_dict.keys())
	genes_from_qtl = set(gene_list)
	
	present_genes = genes_from_file & genes_from_qtl
	missing_genes = genes_from_qtl - present_genes
	
	
	return missing_genes, present_genes
			


	
def annotate_from_csv(data_dict, gene_list):
	"""
	This data comes from a csv file
	data is a dictionary
	each key is a gene name
	the corresponding value is a list of go terms
	
	The output are two dictionaries for genes inside and outside the QTL's
	"""
	
	#Create two new dictionaries; 
	# 1: for all the genes inside the found QTL's
	# 2: for all the genes outside the found QTL's
	gene_inside_qtl_dict = {}
	gene_outside_qtl_dict = {}
	
	#Iterate through the complete data_dict
	#Note that we have to go through ALL the genes in order to
	#populate both dictionaries
	for gene in data_dict:
		
		if gene in gene_list:
			
			#If the gene from the gene list is present in the data_dict
			#then add that gene as key to the new dict with 
			#corresponding value
			#the value is already a list of go annotations that are
			#predicted for that gene
			gene_inside_qtl_dict[gene] = data_dict[gene]
		
		if gene not in gene_list:
			
			#If the gene is not inside the QTL's then it must be outside
			gene_outside_qtl_dict[gene] = data_dict[gene]
	
	#Calculate the totals here so that they can be used to populate
	#The contingency tables very fast
	tot_in_qtl = len(gene_inside_qtl_dict)
	tot_out_qtl = len(gene_outside_qtl_dict)

	
	return gene_inside_qtl_dict, gene_outside_qtl_dict, tot_in_qtl, tot_out_qtl
	

def unique_GO_list(gene_in_qtl_dict):
	"""
	Gene_dict contains a gene name for each key
	The values are lists of go terms for each gene inside the QTL
	
	Extract the go terms from the dictionary values
	And create a list of unique go terms
	"""
	
	golist = []
	
	for gene in gene_in_qtl_dict:
		
		#Basicly for each go term in the value...
		for go in gene_in_qtl_dict[gene]:
			
			#make sure the go term is NOT already in the list
			#because each go term should be unique
			if go not in golist:
				golist.append(go)
				
			else:
				continue
			
	return golist
	
	


########################################################################
########################Prepare contingency table#######################
########################################################################



def lets_see_tables(yesno_list):
	"""
	Just for visualizing the contingency tables...
	"""
	
	c_array = []
	i = 1
	for info in yesno_list:
		#extract go term
		go = info[0]
		#extract contingency table
		
		
		total_go = info[1] + info[2]
		total_nogo = info[3] + info[4]
		total_qtl = info[1] + info[3]
		total_noqtl = info[2] + info[4]
		
		total_total1 = total_qtl + total_noqtl
		total_total2 = total_go + total_nogo
		
		
		c_table = [[i], [info[1], info[2], total_go], [info[3], info[4], total_nogo],[total_qtl, total_noqtl, total_total1, total_total2]]
		
		c_array.append(c_table)
		
		i += 1
		
	print "number of c tables=" , i
		
	return c_array




def make_go_array(gene_in_qtl_dict):
	"""
	Make go_array from dict values
	The values are lists of go terms
	
	Note that individual GO terms can be present in the list multiple 
	times. This is exactly the point
	"""

	go_array = []
	
	for gene in gene_in_qtl_dict:
		
		go_array.append(gene_in_qtl_dict[gene])


	return go_array




def flatten_array(go_array):
	"""
	Flatten one level of nesting so that an array becomes a list
	"""

	flattened_list = []
	for go_list in go_array:
		
		for go in go_list:
			
			flattened_list.append(go)

	return flattened_list
	
	
def count_all_goterms(go_list):
	"""
	Count all go terms in the list in one swift subroutine
	The Counter() function comes from pythons collections module
	
	It counts all values in a list and spits out a dictionary with
	(unique) values from the list as keys and the number of times they
	appear in the list as values
	"""

	counted_go_dict = Counter(go_list)

	
	return counted_go_dict
	



def populate_contingency_table(go_in_qtl_dict, go_out_qtl_dict, golist, total_in_qtl, total_out_qtl):
	"""
	Try to find a faster way to calculate the contingency tables
	-edit: success, went from 10 seconds to 0.01
	
	c_table n: [[0], [1], [2], [3], [4] ]
	
	[0] = 'GO annotation'
	[1] = number of genes with 'GO term' inside all found qtl's
	[2] = number of genes with 'GO term' outside all found qtl's
	[3] = number of genes without 'GO term' inside all found qtl's
	[4] = number of genes without 'GO term' outside all found qtl's
	
	
	c_table layout:
	
			QTL			noQTL
			________________________________	
	GO		|	n11		|	n12		|		|
	noGO	|	n21		|	n22		|		|
	total	|	n31		|	n32		|	n33	|
			|___________|___________|_______|		
	
	"""
	
	c_array = []

	#The marginal totals that are the same for every contingency table
	n31 = total_in_qtl
	n32 = total_out_qtl
	n33 = n31 + n32
	
	for go in golist:
		
		if go in go_in_qtl_dict and go in go_out_qtl_dict:
			
			#Thanks to the known marginal totals we can abuse the degrees
			#of freedom; Once n11 and n12 are known, a simple reduction
			#will tell us the values of n21 and n22
			n11 = go_in_qtl_dict[go]
			n12 = go_out_qtl_dict[go]
			n21 = n31 - n11
			n22 = n32 - n12
			
			c_table = [go, n11, n12, n21, n22]
		
			c_array.append(c_table)
	

	
	return c_array, n33
	


	
###############################################################################
###########################Fisher###Test#######################################
################################Scipy##########################################
###############################################################################



def fish_for_python(data):
	"""
	Using scipy to do the fisher exact test on the created contingency tables
	
	The speed of doing the fishers exact test was compared in scipy and 
	imported R objects. It was found that the test in scipy gave the same results
	but at 10 times the speed.

	The fisher test is performed on each contingeny table
	Each contingency table is a measurement for a single go term.
	
	The output is a go term and the corresponding calculated p-values for
	both alternative hypothesis (Even though we only care about the RPV)
	
	##################################################
	
	Input data structure:
	
			inQTL		outQTL
	GO		info[1]		info[2]
	no-GO	info[3]		info[4]
	
	Example: 
	
	data = [['GO:0000003', 76, 409, 3152, 16392]]
	
	Contingency table:
					inQTL	outQTL	row total
	GO				76		3152	3228
	no-GO			409		16392	16801
	column total 	485		19544	20029
	
	##################################################
	
	The question asked is: 
	
	Knowing that 3228 genes of these 20029 genes are genes with a GO term, 
	
	and that 485 genes of the 20029 genes are inside a QTL,
	 
	and assuming that its equally likely for genes inside and outside 
	the QTL to have a GO term, 
	
	what is the probability that these 3228 genes would be so unevenly 
	distributed between the outside and inside of the QTLs?
	
	IF we were to choose 3228 genes at random, what is the probability that
	3152 genes or more of them would be among the 19544 genes outside the QTL?
	
	and only 76 or fewer from among 485 genes inside the QTL?
	
	"""
	
	i = 1

	fisher_output = []
	
	for info in data:
		#extract go term
		go = info[0]
		#extract contingency table
		c_table = [[info[1], info[2]], [info[3], info[4]]]
		
		#perform two alternative fisher exact test
		fu_oddsratio, fu_pvalue = stats.fisher_exact(c_table, alternative='greater')
		
		#Round the p values to 7 decimals
		fu_pvalue2 = round(fu_pvalue, 7)
		
		individual_fisher_score = [i, go, fu_pvalue2]
		
		fisher_output.append(individual_fisher_score)
		
		
		i += 1
		
		#[i, go, fl_pvalue, fu_pvalue]
		
	
	return fisher_output
	
	
def extract_significant_result(fisher_input, alpha):
	"""
	Only allow the results with a pvalue equal or below alpha to be returned
	Because those are to be considered significant results
	"""
	
	significant_output = []
	
	#info = [i, go, fu_pvalue]
	for info in fisher_input:
		if info[2] <= alpha:
			
		
			significant_output.append(info)
			
	
	return significant_output
	
	
###############################################################################
######################Multiple#Testing#########################################
############################through R##########################################
###############################################################################


def multiple_test(data):
	"""
	Finaly we need R's help to perform some statistics..
	But hopefully not for too long ;)
	
	input = [i, go, fu_p_value[0]]
	
	The p.adjust method is done on the RPV's of the fisher test, 
	since this value is more relevant.
	
	output = [i, go, fu_p_value[0], adjusted(fu_p_value[0]) ]
	
	The multiple testing procedure with the Benjamini/Hochberg method is
	not yet viable in scipy, therefor it is done through imported R functions. 
	
	"""
	adjusted_output = []
	fu_p_value_list = []
	
	#store all the p values derived from the "greater" alternative fisher test
	#in a new list
	for info in data:
		fu_p_value_list.append(info[2])
	
	#transform the list into vector so that R also understands whats
	#going on
	fu_p_value_vector = ro.FloatVector(fu_p_value_list)
	
	#perform the FDR method on the pvalue vector
	cor_data = ro.r["p.adjust"](fu_p_value_vector, method = "BH")
	
	Rround = ro.r['round']
	
	#Iterate through the list of lists and add the q value to each list
	for i in range(0, len(data)):
		
		Rdata = Rround(cor_data[i], digits = 7)
		individual_adjusted = [data[i][0], data[i][1], data[i][2] , Rdata[0]]
		
		adjusted_output.append(individual_adjusted)
		
	return adjusted_output
	
	
	
###############################################################################
############################Post Processing#####################################
####################AaltJan has some demands!###################################
###############################################################################
	




def make_qtl_go_dict(gene_in_qtl_dict, qtl_gene_dict):
	"""
	Make a new dictionary with QTLs as keys.
	replace the values of genelists with corresponding golists
	Resulting in a dictionary of QTLs with a list of the found go terms
	for each QTL
	"""
	
	qtl_go_dict = {}
	
	#Each of these qtls contain a gene_list which will be iterated through
	for qtl in qtl_gene_dict:
		
		gene_loc_list = qtl_gene_dict[qtl]
		temp_go_value = []
		
		#Iterate through gene_list
		for gene_loc in gene_loc_list:
			#for gene in qtl related gene_list
			gene = gene_loc[0]
			
			if gene in gene_in_qtl_dict:
				
				#add golist from the value to the existing golist
				temp_go_value.extend(gene_in_qtl_dict[gene])
			
		
		qtl_go_dict[qtl] = temp_go_value

	
	return qtl_go_dict
				


def post_process_A(enriched_golist, qtl_go_dict, valA):
	"""
	Check for each GO term if it is present in the QTLs
	gene function must be present in more than 50% of the QTL regions of a trait
	"""
	
	approved_golist_A = []
	
	#Iterate through the go terms that were enriched for the genes
	#in the found QTLs
	for go_info in enriched_golist:
		
		go = go_info[1]
		is_go_in_qtl_count = 0.0
		total_qtl = float(len(qtl_go_dict))
		
		#For each QTL in the dictionary, check if a go term is present
		for qtl in qtl_go_dict:
			
			if go in qtl_go_dict[qtl]:
				
				is_go_in_qtl_count += 1.0
		
		go_fraction_score = round(is_go_in_qtl_count / total_qtl, 7)
		
		
		#valA can be given on the website
		#If no value was given, valA will be 0.5
		if go_fraction_score >= valA:
			
			go_info.append(go_fraction_score)
			approved_golist_A.append(go_info)
			
	return approved_golist_A





def post_process_B(approved_golist_A, gene_in_qtl_dict, gene_out_qtl_dict, valB):
	"""
	(B)gene functions must not be predicted for more than 1% of all genes
	
	"""
	
	#Count for how many genes in the QTL a go term from golist is predicted
	#Need approved_golist_A + gene_in_qtl_dict
	
	approved_golist_B = []
	n_genes_in_qtl = len(gene_in_qtl_dict)
	n_genes_out_qtl = len(gene_out_qtl_dict)	
	total_number_genes = float(n_genes_in_qtl + n_genes_out_qtl)
	print "total genes:", total_number_genes
	
	#Inherit the SortedDict which will enable the keys to remain sorted
	go_gene_dict_full = SortedDict()
	go_gene_dict_reduced = SortedDict()
	
	for go_info in approved_golist_A:
		
		go = go_info[1]
		go_count = 0.0
		gene_list_in = []

		
		for gene in gene_in_qtl_dict:
			
			if go in gene_in_qtl_dict[gene]:
				go_count += 1.0
				#only add a gene if its found in a QTL
				gene_list_in.append(gene)				
		
		for gene in gene_out_qtl_dict:
			
			if go in gene_out_qtl_dict[gene]:
				go_count += 1.0				
			
				
		go_fraction_score = round(go_count / total_number_genes, 7)

		#valB can be given on the website
		#If no value was given, valB will be 0.01		
		if go_fraction_score <= valB:

			go_info.append(go_fraction_score)
			go_info_tuple = tuple(go_info)
			go_gene_dict_full[go_info_tuple] = gene_list_in
			
			go_reduced = go_info_tuple[1]
			go_gene_dict_reduced[go_reduced] = gene_list_in
			
	return go_gene_dict_full , go_gene_dict_reduced


			
		
			
###############################################################################
#########################Create url Links######################################
###############################################################################




def GO_link_creator(GOterm):
	"""
	This function will prepare hyperlinks for the identified GO terms
	to be displayed
	"""
	
	#go_list = approved_godict.keys()
	golink_dict = {}
	main_url = "http://www.ebi.ac.uk/QuickGO/GTerm?id="
	
	#for go_info in go_list:
		
		#go = go_info[1]
		#link = main_url + go
		
		#golink_dict[go] = link
		
		
	link = main_url + GOterm
	
	go_link_tuple = (link, GOterm)
		
	return go_link_tuple





def gene_link_creator(gene):
	"""
	Take a trait and create a link to a website
	
	Output is a list with link and gene name
	"""
	
	main_url_begin = "http://www.arabidopsis.org/servlets/Search?type=general&search_action=detail&method=1&show_obsolete=F&name="
	main_url_end = "&sub_type=gene&SEARCH_EXACT=4&SEARCH_CONTAINS=1"

		
	link = main_url_begin + gene + main_url_end
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




def get_genes_with_go_from_qtl(qtl_gene_dict, go_dict):
	"""
	Get the genes from each QTL for which a function was predicted
	
	In order to do this make a set() of the gene lists in go_dict
	And also make a set() of the gene list in qtl_gene_dict
	
	Compare the sets 
	Extend the gene lists
	and the ones that are intersecting should remain
	
	This also has to be done for the gene list per go term
	"""
	
	qtl_gogenes_dict = SortedDict()
	
	for qtl in qtl_gene_dict:
		
		intersecting_genes = []
		gene_list111 = qtl_gene_dict[qtl]
		gene_list11 = [str(i[0]) for i in gene_list111]
		gene_list1 = set(gene_list11)
		
		for go in go_dict:
			
			genelist2 = set(go_dict[go])
			temp_list = gene_list1 & genelist2
			
			intersecting_genes.extend(temp_list)
		
		qtl_gogenes_dict[qtl] = sorted(intersecting_genes)
		
	return qtl_gogenes_dict
		
	


def link_the_linked_links(qtl_gene_dict, go_dict):
	"""
	Prepare the display of the gene list for each go term 
	Divide the genelist of each go term into presence in QTLs
	"""
	
	go_genes_qtl_dict = SortedDict()
	
	for go in go_dict:

		
		#This will be the key in the dictionary
		link_go_tuple = GO_link_creator(go)

		
		intersecting_array = []
		go_gene_set = set(go_dict[go])
		#print sorted(go_gene_set)
		for qtl in qtl_gene_dict:
			
			gene_list111 = qtl_gene_dict[qtl]
			gene_list11 = [str(i[0]) for i in gene_list111]
			qtl_gene_set = set(gene_list11)
			intersecting_genes = go_gene_set & qtl_gene_set
			
			intersecting_gene_link_list = []
			sorted_intersecting_genes = sorted(intersecting_genes)
			
			for gene in sorted_intersecting_genes:
				
				intersecting_gene_link_list.append(gene_link_creator(gene))


				
			intersecting_array.append(sorted(intersecting_gene_link_list))
			
			
		go_genes_qtl_dict[link_go_tuple] = intersecting_array
		
	return go_genes_qtl_dict






###############################################################################
#################################Display#######################################
##################################Output#######################################
###############################################################################

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

		
	#The values 'trait_name' and 'cutoff_name' are given in the template
	if request.GET.get('trait_name') and request.GET.get('cutoff_name') and request.GET.get('experiment_name'):
		if 'gxe' in request.GET:
			gxe_boolean = True
		if 'gxe' not in request.GET:
			gxe_boolean = False

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

		trait_function = get_trait_description(trait_name)
		
		#Create url for GraphView
		Graph_url_base = request.build_absolute_uri(reverse('display_graph'))
		Graph_url_data = "?trait_name=%s&cutoff_name=%f&experiment_name=%s" %(trait_name, cutoff_name, experiment_name)
		Graph_url = Graph_url_base + Graph_url_data
		

		traitquery = Query(trait_name, gxe_boolean, experiment, cutoff_name, fisher_alpha_name, postA_name, postB_name)
		
		traitquery.marker_logp_list()
		traitquery.add_chromosome_and_position()
		traitquery.retrieve_logp_above_cutoff(

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
			filename = 'At_geneGOlists.csv'
			module_dir = os.path.dirname(__file__)  # get current directory
			file_path = os.path.join(module_dir, 'documents', filename)
			
			#read the rows from the csv file mentioned above in filename
			data_dict = read_once_csv(file_path)
			
			absent_genes, present_genes = segregate_gene_list(data_dict, gene_list)


			print "number of genes in gene_list = %d" %len(gene_list)
			print "number of genes in data_dict =", (len(data_dict) - len(gene_list))	
			
			gene_inside_qtl_dict, gene_outside_qtl_dict, tot_in_qtl, tot_out_qtl = annotate_from_csv(data_dict, present_genes)
			print "missings = %d" % len(absent_genes)
			
			print "number of genes in qtl = %d and %d" %(len(gene_inside_qtl_dict), len(present_genes))
			print "number of genes out qtl =", len(gene_outside_qtl_dict)
		
			golist_unique = unique_GO_list(gene_inside_qtl_dict)	
			

			##############################################################
			##############################################################		


			#Count all genes with and without GO annotations
			#Count all genes inside and outside the QTLs
			
			#Prepare for counted genes inside qtl:
			golist_array_in = make_go_array(gene_inside_qtl_dict)
			golist_in = flatten_array(golist_array_in)
			c_go_in_qtl_dict = count_all_goterms(golist_in)
			
			
			#Prepare for counted genes outside qtl:
			golist_array_out = make_go_array(gene_outside_qtl_dict)
			golist_out = flatten_array(golist_array_out)
			c_go_out_qtl_dict = count_all_goterms(golist_out)
			
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
			enriched_golist = multiple_test(significant_info)
			
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
			godict_full, godict_reduced = post_process_B(approved_golistA, gene_inside_qtl_dict, gene_outside_qtl_dict, postB_name)
			toc = time.clock()
			print "Time to approve enrichment B is %f seconds" % (toc-tic)
			print "length of godict = %d" %len(godict_full)
			
			
			
			#GO links:
			link_for_trait = gene_link_creator(trait_name)
			
			
			#go_gene_link_dict = link_the_linked_links2(godict_reduced)
			
			qtl_gogenes_dict = get_genes_with_go_from_qtl(gene_dict, godict_reduced)
			
			go_genes_qtl_dict = link_the_linked_links(gene_dict, godict_reduced)
			
			
			
			#c_array_display = lets_see_tables(c_array)
			
			

			

			
					
			stats_dict = {
						"link_for_trait" : link_for_trait,
						"Graph_url" : Graph_url,
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
		trait_list = trait_query.split()
		
		combined_stats = []
		
		#translate cuttoff
		cutoff_name_D = request.GET.get('cutoff_name')
		cutoff_name = float(cutoff_name_D)
		
		experiment_name = request.GET.get('experiment_name').strip()
		
		
		#get list of all traits
		#traits = Gene.objects.values_list('locus_identifier', flat=True).order_by('locus_identifier')

		#find number of qtl for trait
		for trait in trait_list:
			
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
		response = HttpResponse(imgData.getvalue(), mimetype = 'image/png')
		
		#response = HttpResponse(mimetype = 'image/png')
		#canvas.print_png(response)
		
		toc = time.clock()

		print "time for graph construction is %f" %(toc-tic)
		return response
			

	else:
		exps = Experiment.objects.all().values_list('experiment_name',flat=True)
		return render_to_response('wouter/search_graph_trait.html',{'exps':exps})
		#return render_to_response("wouter/search_graph_trait.html")
	




def StoreCsvView(request):
	"""
	
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
	if request.GET.get('trait_name') and request.GET.get('cutoff_name'):
		
	
		trait_name_u = request.GET.get('trait_name').strip()
		trait_name = trait_name_u.encode("utf-8")
	
		cutoff_name_D = request.GET.get('cutoff_name').strip()
		cutoff_name = float(cutoff_name_D)
	
		trait_function = get_trait_description(trait_name)

		l1 = marker_logp_list(trait_name, gxe_boolean)
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
		return render_to_response('wouter/search_individual_trait.html')
		
		
		
	
