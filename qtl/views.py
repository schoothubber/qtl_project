import os
import time
import csv

import rpy2.robjects as ro
from rpy2.robjects.vectors import FloatVector
import numpy
import scipy
from scipy import stats

from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseRedirect
from django.http import HttpRequest as request
from django.template import RequestContext

from .models import Experiment,Gene,Marker,LOD,Parent,RIL,Metabolite,MParent,MRIL,MLOD



###############################################################################
#################################WOUTERS#######################################
###################################CODE########################################
###############################################################################


def marker_logp_list(trait, gxe_boolean):
	'''
	Create a list of tuples containing markers with corresponding -log p value.
	Both are subscribed to the requested trait.
	The markers and -log p values are ordered based on the markers' location on the
	chromosomes and their genetic mapping in centiMorgan.
	 
	The list is enumerated, meaning that an integer is added to
	the tuple. The integer ranges from 0 to len(pair_list)-1
	The enumeration will eventually help to identify marker regions
	in the check_if_markers_consecutive function
	'''
	
	#First get all the marker names from the database and order them
	#according chromosomal position and genetic mapping
	markers = Marker.objects.values_list('marker_name',flat=True).order_by('marker_chromosome','marker_cm')
	
	marker_logp_tuple = ()
	tuple_list = []
	 
	# iterate of the list of markers and get the -log p value for each marker
	for i in range(0, len(markers)):
		
		try:
			#Because there is only one object that will match the query
			#the get() method is used on a Manager (in this case objects)
			#to return the object directly
			#gxe differentiates between normal gene expression (False)
			#and environmental permutation (True)
			logp1 = LOD.objects.get(locus_identifier = trait, marker_name = markers[i], gxe = gxe_boolean)
			
			#From the returned object extract the value from the attribute LOD_score
			logp2 = logp1.LOD_score
			
			#Change the type from "Class decimal.Decimal" to float
			logp3 = float(logp2)#the actual -log p value
			
			#add the i for enumeration
			#perform encode("utf-8") on each marker to reduce the sice
			marker_logp_tuple = (i, markers[i].encode("utf-8"), logp3)
			tuple_list.append(marker_logp_tuple)
	        
		except LOD.DoesNotExist:
			pass
			#logp3 = None


	#checked_tuple_list = []
	#for tup in tuple_list:
		
		#if tup[2] equals to None then the gene name does not exist in the database
		#there for it cannot be used fo further manipulation
		#if tup[2] != None:
			#checked_tuple_list.append(tup)
	
	
	return tuple_list
		#trait, [(int, marker name, -logp)]
		

	
	
###############################################################################
###############################################################################

		
def add_chromosome_and_position(tuple_list): #[(int, marker, -logp)]
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


###############################################################################
###############################################################################


def retrieve_logp_above_cutoff(tuple_list, cutoff):
	'''
	Take all the Markers and -logp values above a certain value,
	(Usually a cutoff of 2.3 is seen as significant)
	'''

	cutoff_list = []
	# Iterate through the tuple_list and exclude all tuples with a
	#-logp value lower than the cutoff
	for info_tuple in tuple_list:#for each tuple (int, marker, lodscore)
		
		
		#Use abs() because there are negative -logp values present
		#The minus signs were added to indicate the differential
		#expression from the other allel
		# + indicates Bay and - indicates Sha
		if abs(info_tuple[2]) >= cutoff:
			
				
			cutoff_list.append(info_tuple)
		else:
			continue
			
	return cutoff_list #(int, marker, -logp)
	
	

###############################################################################
###############################################################################


def check_region(cutoff_list):
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
   
	
###############################################################################
###############################################################################


def iterate_marker_tuple(marker_region_list):
	"""
   Iterate over a list containing lists with either one or more
   tuples. 
   A single tuple will be identified as the maximum in the group.
   A list of 2 or more tuples will be fed to the function 
   highest_tuple. 
   """
    
	maximum_markers = []
    
    #Iterate over the lists containing the tuple(s)
	for region in marker_region_list:
	
		maximum_markers.append(highest_tuple(region))
	
	return maximum_markers


###############################################################################

	
def highest_tuple(region):
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



###############################################################################
###############################################################################
        
        
def set_region_from_adjacent_markers(maximum_list, marker_list):
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
			chr_max = get_chromosome_max(maximum_list[i])
			
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
			corr_start, corr_end = check_regions_on_chromosome(region_start, region_end, marker_list)
			
			actual_region_list.append([corr_start, corr_end])
	
			
    return actual_region_list


###############################################################################


def get_chromosome_max(marker_tuple):
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
		
		
###############################################################################


def check_regions_on_chromosome(region_start, region_end, marker_list):
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
			chr_max = get_chromosome_max(region_start)
			corr_end = (None, None, None, region_start[3], chr_max)
			return corr_start, corr_end
		
		#if the topmarkers chromosome number is equal to the end marker
		#then keep the end marker as it is and alter the start marker
		elif revert_tuple[3] == region_end[3]:
			#The region start should be on the same chromosome as the
			#region end.
			#The None values should not matter, this new tuple only has
			#a physical position on a chromosome
			corr_start = (None, None, None, region_end[3], 0)
			corr_end = region_end
			return corr_start, corr_end
			
	
###############################################################################
###############################################################################


def combine_marker_region_data(trait, region_tuple_list):
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
		combined_tuple = (trait, region[0][3], region[0][4], region[1][4])
		combined_data.append(combined_tuple)
		
		
	#combined_data(trait, chromosome , start physical location, end physical location)
	return combined_data
		

###############################################################################
###############################################################################


def find_genes_in_regions(region_list):
	"""
	The region_list will contain tuples with the following information:
	(trait, chromosome number, marker start position, marker end position)
	
	This information is used to collect lists of genes that are positioned
	on each region.
	[gene name, gene start, gene end]
	The tuple and gene lists will be placed together in a dictionary with
	the tuple as key and the list as value
	"""
	
	gene_dict = {}
	
	complete_gene_list = Gene.objects.values_list("locus_identifier","chromosome", "start", "end")
	
	for region in region_list:
		#Only look at the genes on a specific chromosome
		chromosome_spec_list = complete_gene_list.filter(chromosome = region[1])
		#Exclude all genes with basepair locations below the first marker
		remove_front_genes = chromosome_spec_list.exclude(start__lt = region[2])
		#Exclude all genes with basepair locations above the last marker
		remove_end_genes = remove_front_genes.exclude(end__gt = region[3])
		
		#normalize the database data
		normal_list = normalize_object_data(remove_end_genes)
		
		#add each region and list of genes to the dictionary gene_dict
		#with the tuple as a key and the gene list as a value
		gene_dict[region] = normal_list
		
	return gene_dict
	
	
###############################################################################


def normalize_object_data(data_list):
	"""
	The extracted data from the database needs to be parsed 
	data_list has the structure 
	[gene name,chromosome number, gene start, gene end]
	With data types that were derived from the model classes
	The output should become 
	[gene name, gene start, gene end]
	With normalized data types
	"""	
	normal_list = []
	
	for data in data_list:
		normal_list.append([data[0].encode("utf-8"), int(data[2]), int(data[3])])
	
	return normal_list


###############################################################################
###############################################################################


def read_distinct_genes(genedict):
	"""
	This function gives all found genes from all regions an annotation
	using the GO annotations from a file.
	All genes from all regions will be placed in one big list.
	The ..
	"""
	
	#read all the genes from the found regions into a new list
	#This list contains only genes which have to be unique
	#i.e. no gene can be mentioned more than once
	gene_list = []
	
	for key in genedict:
		
		for genedata in genedict[key]:
			
			if genedata[0] not in gene_list:
				gene_list.append(genedata[0])
			
			
	return gene_list

###############################################################################
############################Implement##########################################
########################AaltJan's##scripts#####################################
###############################################################################
#############################-EDIT#############################################
###################AaltJan's#scripts###has been################################
###########################################replaced############################
###############################################################################



def read_once_csv(filename):
	"""
	Read a csv file and place the contents in a dictionary
	
	The csv file contains rows of data
	In each row the first index holds a gene name
	The rest of the row holds the GO annotations
	"""
	
	tic = time.clock()
	
	#data = genfromtxt(filename, delimiter = ',')
	data_dict = {}
	
	with open(filename, 'rb') as csvobject:
		reader = csv.reader(csvobject)
		for row in reader:
			
			#key is the gene name
			key = row[0]
			#value is a list of go terms
			value = row[1:]
			
			data_dict[key] = value
		
	
	toc = time.clock()
	stopwatch = toc-tic
	
	return data_dict, stopwatch



	
def annotate_from_csv(data_dict, gene_list):
	"""
	This data comes from a csv file
	data is a dictionary
	each key is a gene name
	the corresponding value is a list of go terms
	"""
	tic = time.clock()
	
	#Create two new dictionaries; 
	# one for all the genes inside the found QTL's
	# other one for all the genes outside the found QTL's
	gene_inside_qtl_dict = {}
	gene_outside_qtl_dict = {}
	
	#Iterate through the complete gene list
	for key in data_dict:
		
		if key in gene_list:
			
			#If the gene from the gene list is present in the data_dict
			#then add that gene as key to the new dict with 
			#corresponding value
			#the value is already a list of go annotations that are
			#linked to that gene
			gene_inside_qtl_dict[key] = data_dict[key]
		
		if key not in gene_list:
			
			#If the gene is not inside the QTL's then it must be outside
			gene_outside_qtl_dict[key] = data_dict[key]

	how_many_in = len(gene_inside_qtl_dict)
	how_many_out = len(gene_outside_qtl_dict)

	toc = time.clock()
	stopwatch = toc-tic
	
	return gene_inside_qtl_dict, gene_outside_qtl_dict, how_many_in, how_many_out, stopwatch
	

def unique_GO_list(gene_dict):
	"""
	Gene_dict contains a gene name for each key
	The values are lists of go terms for each gene
	
	Extract the go terms from the dictionary values
	And create a list of unique go terms
	"""
	
	golist = []
	
	for key in gene_dict:
		
		#Basicly for each go term in the value...
		for go in gene_dict[key]:
			
			#make sure the go term is NOT already in the list
			#because each go term should be unique
			if go not in golist:
				golist.append(go)
				
			else:
				continue
			
	return golist
	
	


#def get_geneGO_outside_qtl(data_dict, qtl_gene_dict):
	#"""
	#What goes in is:
	#data_dict[gene_name] = [list of go terms] (of the whole genome)
	#qtl_gene_dict[gene_name] = [list of go terms] (inside the qtls)
	
	#What comes out is:
	#outsideqtl_dict[gene_name] = [list of go terms] (outside the qtls)
	

	#"""


	#outsideqtl_dict = data_dict
	
	##Check each gene from the ca. 20.000 rows
	#for key in qtl_gene_dict:
		
		##If this gene was not found inside the QTL's
		##Then add it to the list of genes outside the QTL's
		#if key in outsideqtl_dict.keys():
			
			#del outsideqtl_dict[key]
			
			
	#return outsideqtl_dict
		
	


########################################################################
########################Prepare contingency table#######################
########################################################################

def make_yesno_list(golist, seldict, alldict):
	"""
	The code in this function is based on script written by Aalt-Jan van Dijk
	
	Output yesno_list: [[0], [1], [2], [3], [4] ]
	
	[0] = 'GO annotation'
	[1] = number of genes with 'GO term' inside all found qtl's
	[2] = number of genes without 'GO term' inside all found qtl's
	[3] = number of genes with 'GO term' outside all found qtl's
	[4] = number of genes without 'GO term' outside all found qtl's
	
	seldict contains all genes, with corresponding GO terms, inside the qtl's
	whereas
	alldict contains all genes, with corresponding GO terms, outside the qtl's
	"""
	
	yesno_list = []
	
	#Iterate through each GO term in the golist
	for go in golist: # [0]
		
		in_yes = 0 # [1]
		in_no = 0 # [2]
		
		out_yes = 0 # [3]
		out_no = 0 # [4]
		
		#Iterate through each gene from seldict
		#calculate output for [1] and [2] 
		for gene in seldict.keys():
			
			if go in seldict[gene]:
				in_yes += 1
			else:
				in_no += 1
				
		
		#Iterate through each gene from alldict
		#calculate output for [3] and [4] 
		for gene in alldict.keys():
			
			if go in alldict[gene]:
				out_yes += 1
			else:
				out_no += 1
				
		go_yesno = [go, in_yes, in_no, out_yes, out_no]
		yesno_list.append(go_yesno)
		
	return yesno_list



###############################################################################
###########################Fisher###Test#######################################
#############################   R   ###########################################
###############################################################################


	
def fisher_test(yesno_list):
	"""
	Takes the output of make_yesno_list. For the datastructure of input
	see DOC string of make_yesno_list
	
	The data is fed to robjects using the rpy2 package
	The fisher test is performed for each go term.
	
	The output is a go term and the corresponding calculated p-values for
	both alternative hypothesis
	
	##################################################
	
	Input data structure:
	
			inQTL		outQTL
	GO		info[1]		info[3]
	no-GO	info[2]		info[4]
	
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
	tic = time.clock()
	
	fisher_output = []
	
	#Create a python function called 'go_fish' to call an R function
	#to perform the fisher test
	go_fish = ro.r['fisher.test']

	i = 1
	
	for info in yesno_list:
		#stash the go term here for later use
		go = info[0]
		
		#Create a vector with the values from info
		a = ro.IntVector([info[1], info[3], info[2] ,info[4]])
		#Turn it into a matrix (2x2) for the fisher test
		aa = ro.r['matrix'](a, nrow = 2)
		
		#Perform the fisher test
		# "alternative" indicates the alternative hypothesis and must be 
		#one of "two.sided", "greater" or "less"
		fl = go_fish(aa,alternative = "l")
		fu = go_fish(aa,alternative = "g")

		#Get the p-values from both fisher tests
		#The fisher output is an R vector-like object
		#Items can be accessed with:
		#the delegators rx or rx2
		Rround = ro.r['round']
		fl_p_value = Rround(fl.rx2('p.value'), digits = 7)
		fu_p_value = Rround(fu.rx2('p.value'), digits = 7)
		
		#Place the p-values for each go term in a list:
		#[i, go, fl:p-value, fu:p-value]
		individual_fisher_score = [i, go, fl_p_value[0], fu_p_value[0]]
		
		fisher_output.append(individual_fisher_score)
		
		i += 1
	
	toc = time.clock()
	
	stopwatch = toc-tic
		
	return fisher_output
	
	
	
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
	"""
	
	i = 1

	fisher_output = []
	
	for info in data:
		#extract go term
		go = info[0]
		#extract contingency table
		c_table = [[info[1], info[3]], [info[2], info[4]]]
		
		#perform two alternative fisher exact test
		fl_oddsratio, fl_pvalue = stats.fisher_exact(c_table, alternative='less')
		fu_oddsratio, fu_pvalue = stats.fisher_exact(c_table, alternative='greater')
		
		#Round the p values to 7 decimals
		fl_pvalue2 = round(fl_pvalue, 7)
		fu_pvalue2 = round(fu_pvalue, 7)
		
		individual_fisher_score = [i, go, fl_pvalue2, fu_pvalue2]
		
		fisher_output.append(individual_fisher_score)
		
		
		i += 1
		
		#[i, go, fl_pvalue, fu_pvalue]
		
	
	return fisher_output
	
	
def extract_significant_result(fisher_input, alpha):
	"""
	Only allow the results with a pvalue equal or below alpha to be returned
	"""
	
	significant_output = []
	
	#info = [i, go, fl_pvalue, fu_pvalue]
	for info in fisher_input:
		
		if info[3] <= alpha:
			
		
			significant_output.append(info)
			
	
	return significant_output
	
	
###############################################################################
######################Multiple#Testing#########################################
############################through R##########################################
###############################################################################


def multiple_test(data):
	"""
	input = [i, go, fl_p_value[0], fu_p_value[0]]
	
	The p.adjust method is done on the pvalue of the fisher test with the
	greater than alternative, since this value is more relevant.
	
	output = [i, go, fl_p_value[0], fu_p_value[0], adjusted(fu_p_value[0]) ]
	
	The multiple testing procedure with the Benjamini/Hochberg method is
	not yet viable in scipy, therefor it is done through imported R functions. 
	
	"""
	adjusted_output = []
	fu_p_value_list = []
	
	#store all the p values derived from the "greater" alternative fisher test
	#in a new list
	for info in data:
		fu_p_value_list.append(info[3])
	
	#transform the list into vector so that R also understands whats
	#going on
	fu_p_value_vector = ro.FloatVector(fu_p_value_list)
	
	#perform the FDR method on the pvalue vector
	cor_data = ro.r["p.adjust"](fu_p_value_vector, method = "BH")
	
	Rround = ro.r['round']
	
	#Iterate through the list of lists and add the q value to each list
	for i in range(0, len(data)):
		
		Rdata = Rround(cor_data[i], digits = 7)
		individual_adjusted = [data[i][0], data[i][1], data[i][2] ,data[i][3], Rdata[0]]
		
		adjusted_output.append(individual_adjusted)
		
	return adjusted_output



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
		
		tic_t = time.clock() # starting point of total process
		tic = time.clock() # starting point of genelist retrieval
	
		trait_name_u = request.GET.get('trait_name')
		trait_name = trait_name_u.encode("utf-8")
	
		cutoff_name_D = request.GET.get('cutoff_name')
		cutoff_name = float(cutoff_name_D)
	
		l1 = marker_logp_list(trait_name, gxe_boolean)
		l2 = add_chromosome_and_position(l1)
		l3 = retrieve_logp_above_cutoff(l2, cutoff_name)
		
		#If no markers are found above the chosen cutoff, than display this message
		if l3 == None:
			message = "There are no markers with a LOD score higher than %f for trait %s" % (cutoff_name, trait_name)
			bottle = {"message" : message}
			return render_to_response("wouter/no_go.html", bottle)
		
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
			l8 = find_genes_in_regions(l7)
			
			gene_list = read_distinct_genes(l8)
			
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
	
			
			tic = time.clock()
			
			#get the location of the GO annotation file
			filename = 'At_geneGOlists.csv'
			module_dir = os.path.dirname(__file__)  # get current directory
			file_path = os.path.join(module_dir, 'documents', filename)
			
			#read the rows from the csv file mentioned above in filename
			data_dict, timed_read = read_once_csv(file_path)
			
			#Create a dictionary containing all genes inside the found qtls
			#with corresponding GO annotations	
			########qtl_gene_dict, n_annotations, timed_annotate = annotate_from_csv(data_dict, gene_list)
			
			gene_inside_qtl_dict, gene_outside_qtl_dict, how_many_in, how_many_out, stopwatch = annotate_from_csv(data_dict, gene_list)
			
			golist = unique_GO_list(gene_inside_qtl_dict)
			
			#Create a dictionary containing all genes outside the found qtls
			#with corresponding GO annotations
			#################outside_qtl_dict = get_geneGO_outside_qtl(data_dict, qtl_gene_dict)
			
			#Create a contingency table for the Fishers exact test
			#Count all genes with and without GO annotations
			#Count all genes inside and outside the QTLs
			yesno_list = make_yesno_list(golist, gene_inside_qtl_dict, gene_outside_qtl_dict)

			
			
			n_qtl = len(l8) #number of qtl's 
			n_genes = len(gene_list) #total number of genes from all qtl's
			#n_GO = len(l10) #number of GO annotations found from GO list
						
						
			gene_list_dict = {
						#"n_GO" : n_GO,
						"n_qtl": n_qtl,
						"n_genes": n_genes,
						#"l10": l10,
						"l8": l8,
						"yesno_list": yesno_list
						}
						
			toc = time.clock()
			print "Time to annotate genes is %f seconds" % (toc-tic)
				
				
				
			##############################################################
			####################Fishing with SciPy########################
			##############################################################
			
			
			
			#perform the fisher exact test on all created contingency
			#tables
			#fisher_python = fish_for_python(yesno_list)
			fisher_R = fisher_test(yesno_list)
			
			#only allow the results with p values below alpha to pass
			alpha = 0.05
			significant_info = extract_significant_result(fisher_R, alpha)
			
			#perform a multiple test
			add_adjusted_p = multiple_test(significant_info)
			
			
						
			stats_dict = {
						"trait_name" : trait_name,
						"cutoff_name" : cutoff_name,
						#"fisher_python" : fisher_python,
						#"significant_info" : significant_info,
						"add_adjusted_p" : add_adjusted_p
						}
	
			toc_t = time.clock()
			print "Time to perform entire process is %f seconds" % (toc_t-tic_t)
			
			
			#return display_meta(request)
			return render_to_response("wouter/get_go.html", stats_dict)

		
	else:
		return render_to_response('wouter/search_trait.html')







def CheckTraitView(request):
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
		
		
		#get list of all traits
		#traits = Gene.objects.values_list('locus_identifier', flat=True).order_by('locus_identifier')

		#find number of qtl for trait
		for trait in trait_list:
			
			trait_name = trait.encode("utf-8")
			
			#First search for the amount of QTL's
			l1 = marker_logp_list(trait_name, gxe_boolean)
			l2 = add_chromosome_and_position(l1)
			l3 = retrieve_logp_above_cutoff(l2, cutoff_name)
			l4 = check_region(l3)
			l5 = iterate_marker_tuple(l4)
			l6 = set_region_from_adjacent_markers(l5, l2)
			l7 = combine_marker_region_data(trait_name, l6)
			
			#Howmany qtl for this trait?....
			number_of_qtl = len(l7)
			
			#Second search for the amount of GO enrichments
			l8 = find_genes_in_regions(l7)
			gene_list = read_distinct_genes(l8)
			
			number_of_genes = len(gene_list)
			
			
			combined_stats.append([trait_name, number_of_qtl, number_of_genes])
			
			
			print "trait %s has %d qtls with %d genes" % (trait_name, number_of_qtl, number_of_genes)

		toc = time.clock()
		print "Time to perform trait check is %f seconds" % (toc-tic)
	
		stats_dict = {
					"combined_stats" : combined_stats
					}
	
	
		return render_to_response("wouter/get_trait_stats.html", stats_dict)
	
	
	else:
		return render_to_response("wouter/check_trait.html")







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
		

def SearchTraitView2(request):
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
		
		tic_t = time.clock() # starting point of total process
		tic = time.clock() # starting point of genelist retrieval
	
		trait_name_u = request.GET.get('trait_name')
		trait_name = trait_name_u.encode("utf-8")
	
		cutoff_name_D = request.GET.get('cutoff_name')
		cutoff_name = float(cutoff_name_D)
	
		l1 = marker_logp_list(trait_name, gxe_boolean)
		l2 = add_chromosome_and_position(l1)
		l3 = retrieve_logp_above_cutoff(l2, cutoff_name)
		
		#If no markers are found above the chosen cutoff, than display this message
		if l3 == None:
			message = "There are no markers with a LOD score higher than %f for trait %s" % (cutoff_name, trait_name)
			bottle = {"message" : message}
			return render_to_response("wouter/no_go.html", bottle)
		
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
			l8 = find_genes_in_regions(l7)
			
			l9 = read_distinct_genes(l8)
			
			#Place the output into a dictionary
			#Give each key the same name as the variable for convenience
			#The dictionary will be passed to the template
			#Inside the template the key's are used to display the
			#corresponding values
			
			toc = time.clock()
			print "Time to retrieve a gene list is %f seconds" % (toc-tic)
			
			
			
			gene_dict = {
						"trait_name" : trait_name,
						"cutoff_name" : cutoff_name,
						"l2" : l2,
						"l3" : l3,
						"l4" : l4,
						"l5" : l5,
						"l6" : l6,
						"l8" : l8,
						"l9" : l9,
						}



			return render_to_response("wouter/gene_list.html", gene_dict)

		
	else:
		return render_to_response('wouter/search_trait.html')
