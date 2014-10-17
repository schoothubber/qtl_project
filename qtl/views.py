import string
import os
import time

from django.shortcuts import render,get_object_or_404,render_to_response
from django.http import HttpResponse, HttpResponseRedirect
from django.http import HttpRequest as request
from django.template import RequestContext
from django.core.context_processors import csrf

from .models import Experiment,Gene,Marker,LOD,Parent,RIL,Metabolite,MParent,MRIL,MLOD

###############################################################################
#################################WOUTERS#######################################
###################################CODE########################################
###############################################################################

def marker_logp_list(trait):
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
		
		#Because there is only one object that will match the query
		#the get() method is used on a Manager (in this case objects)
		#to return the object directly
		#gxe differentiates between normal gene expression (False)
        #and environmental permutation (True)
        logp1 = LOD.objects.get(locus_identifier = trait, marker_name = markers[i], gxe = False)
        
        #From the returned object extract the value from the attribute LOD_score
        logp2 = logp1.LOD_score
        
        #Change the type from "Class decimal.Decimal" to float
        logp3 = float(logp2)#the actual -log p value
        
        #add the i for enumeration
        #perform encode("utf-8") on each marker to reduce the sice
        marker_logp_tuple = (i, markers[i].encode("utf-8"), logp3)
        tuple_list.append(marker_logp_tuple)
     
    if tuple_list:
		return tuple_list
		#trait, [(int, marker name, -logp)]
		
    else:
		return "There are no markers with -log p scores for this gene"
	
	
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
		#The minus signs were added to indicate down regulation
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
			#The negative signs weer added after the calculation of the 
			#-logp value to indicate down regulation
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
		elif maximum_list[i][0] == len(marker_list):
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



def read_once(filename):
	"""
	
	"""
	
	fileobject = open(filename, 'r')
	data = fileobject.readlines()
	fileobject.close()
	
	return data



def read_annotation(filename, gene_list):
	"""
	The code in this function is based on code written by Aalt-Jan van Dijk
	
	This funtions takes a list containing the total number of genes from
	all found qtl's, plus a list with similar genes and GO annotations.
	This GO list is a file called 'QTL/qtl/documents/pred.out'
	
	The output is a new list containing the genes that were present in 
	the first gene list AND in the GO list, plus the corresponding
	GO annotation
	
	The tic-toc was added to ascertain how long a search takes
	"""

	#read the data from the file containing the following datastructure:
	#[genename, GOterm, some obscure value]
	
	fileobject = open(filename, 'r')
	data = fileobject.readlines()
	fileobject.close()
	
	Gene_GO_list = []
	
	for gene_GO in data:

		s = string.split(gene_GO)
	
		#if the number of gene names in gene_list is more than zero
		if gene_list.count(s[0])>0:
		
			#only store the [genename, GOterm] part of the list
			info = s[:-1]
			Gene_GO_list.append(info)

	return Gene_GO_list


def readsel(gene_GO_list):
	"""
	The code in this function is based on script written by Aalt-Jan van Dijk
	
	gene_GO_list contains all genes with corresponding GO annotations that
	were found in all qtls
	
	seldict will contain each gene from inside the found qtl as key
	and as value a list of GO terms that subscribe to one gene
	
	"""
	#input is output van getpredfun:
	#[genename, Goterm]	
	
	seldict={}
	golist=[]
	
	for gene_GO in gene_GO_list:
	
		gene = gene_GO[0] #gene name
		go = gene_GO[1]   #Goterm
	
		#If the number of Goterms in the golist is zero
		##################if golist.count(go) == 0:
		if go not in golist:
			
			#add the Goterm to the golist
			#so that each Goterm will be unique in this golist
			golist.append(go)
	
		#if the dictionary already has a key named "genename"
		if seldict.has_key(gene):
			
			#tmp becomes the value beloning to seldict with key "gene"
			#the value is a list of Go terms
			tmp = seldict[gene]
	  
		else:
			#make new list
			tmp=[]
			
		# add Goterm to the new list
		tmp.append(go)
	
		#add the new golist as value to a dictionary 
		#with gene name as key 
		seldict[gene] = tmp
	
		# return a dictionary with genes as keys and as value a list of 
		# the Goterms for each gene
		# plus a list with all used Goterms
		
	return seldict, golist
	

def readall(filename, seldict):
	"""
	The code in this function is based on script written by Aalt-Jan van Dijk
	
	Excludelist is a list of the keys from seldict.readsel which contain
	the genes that were found inside the qtls.
	
	filename = ????
	"""
	#Make a new list and populate it with the key's from seldict
	excludelist = []
	for key in seldict:
		excludelist.append(seldict[key])
			
	#Open a file and read its contents
	
	fileobject = open(filename, 'r')
	data = fileobject.readlines()
	fileobject.close()
	
	alldict = {}
	

	#Iterate through each line of text in data
	for i in data:
		s = string.split(i)
		
		gene = s[0]
		go = s[1]
		
		#If this list does not contain this gene
		####################if excludelist.count(gene) == 0:
		if gene not in excludelist:
	  
			# If this dictionary has a key named "gene name"  
			if alldict.has_key(gene):
				
				#tmp becomes the value beloning to alldict with key "gene"
				#the value is a list of Go terms
				tmp = alldict[gene]
    
			else:
				#make a new list named tmp
				tmp = []
				
			# add the go term to the new list
			tmp.append(go)
   
			# set tmp as value with the gene name as key in alldict
			#tmp is a list of go terms
			alldict[gene] = tmp
   
	return alldict

#seldict,golist=readsel(sys.argv[1])
#alldict=readall(sys.argv[2],seldict.keys())

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
			
			if seldict[gene].count(go) > 0:
				in_yes += 1
			else:
				in_no += 1
				
		
		#Iterate through each gene from alldict
		#calculate output for [3] and [4] 
		for gene in alldict.keys():
			
			if alldict[gene].count(go) > 0:
				out_yes += 1
			else:
				out_no += 1
				
		go_yesno = [go, in_yes, in_no, out_yes, out_no]
		yesno_list.append(go_yesno)
		
	return yesno_list
	
def check_yesno_totals(yesno_list):
	"""
	yesno[1] + yesno[2]: Add up the number of genes inside the qtls
	yesno[3] + yesno[4]: Add up the number of genes outside the qtls
	All total1's should be the same
	All total2's should be the same
	"""

	check = []
	
	for yesno in yesno_list:
		total1 = 0
		total2 = 0
		total1 = yesno[1] + yesno[2]
		total2 = yesno[3] + yesno[4]
		check.append([total1, total2])
	
	return check
		

###############################################################################
#################################Display#######################################
##################################Output#######################################
###############################################################################

	
def SearchTraitView(request):
	"""
	When a trait and a cutoff value are given on the website, these values
	are retrieved from the database and fed to this view. 
	They are subsequently manipulated by other functions
	The output is sent to a nother template galled genelist.html
	"""
	#if request.GET.get('gxe') == True:
		#print "gxe = True"
		
	#if request.GET.get('gxe') == False:
		#print "gxe = False"
		
	if request.GET.get('trait_name') and request.GET.get('cutoff_name'):
		
		tic = time.clock()
		
		trait_name_u = request.GET.get('trait_name')
		trait_name = trait_name_u.encode("utf-8")
		
		cutoff_name_D = request.GET.get('cutoff_name')
		cutoff_name = float(cutoff_name_D)
	
		l1 = marker_logp_list(trait_name)
		l2 = add_chromosome_and_position(l1)
		l3 = retrieve_logp_above_cutoff(l2, cutoff_name)
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
		
		###############################
		###############AJ##############
		###############################

		#get the location of the GO annotation file
		filename = 'pred.out'
		module_dir = os.path.dirname(__file__)  # get current directory
		file_path = os.path.join(module_dir, 'documents', filename)
		
		data = read_once(file_path)
				
		l10 = read_annotation(file_path, l9)
		
		seldict, golist = readsel(l10)
		
		alldict = readall(file_path, seldict)
		
		yesno_list = make_yesno_list(golist, seldict, alldict)
		
		check = check_yesno_totals(yesno_list)
		
		
		
		
		
		#Place the output into a dictionary
		#Give each key the same name as the variable for convenience
		#The dictionary will be passed to the template
		#Inside the template the key's are used to display the
		#corresponding values
		gene_dict = {
					"trait_name" : trait_name,
					"cutoff_name" : cutoff_name,
					"l2" : l2,
					"l3" : l3,
					"l4" : l4,
					"l5" : l5,
					"l6" : l6,
					"l8" : l8,
					"l9" : l9
					}
		
		n_qtl = len(l8) #number of qtl's 
		n_genes = len(l9) #total number of genes from all qtl's
		n_GO = len(l10) #number of GO annotations found from GO list
		
		gene_list_dict = {
					"n_GO" : n_GO,
					"n_qtl": n_qtl,
					"n_genes": n_genes,
					"l10": l10,
					"l8": l8,
					"yesno_list": yesno_list,
					"check": check
					}


		toc = time.clock()
		print "in %f seconds" % (toc-tic)

		return render_to_response("wouter/get_go.html", gene_list_dict)
		#return render_to_response("wouter/gene_list.html", gene_dict)

	else:
		return render_to_response('wouter/search_trait.html')
		


		
		

