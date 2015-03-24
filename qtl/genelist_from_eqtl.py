import sys
import os
from sortedcontainers import SortedDict

from qtl.models import Gene,Marker,LOD,Experiment,ExperimentMarker


##Set the path so that the functions in this file can be imported by views.py
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
#sys.path.append('/mnt/geninf15/prog/www/django/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')


def marker_logp_list(trait, gxe_boolean, experiment):
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

	markers_full = ExperimentMarker.objects.filter(experiment_name = experiment)
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
			logp1 = LOD.objects.get(locus_identifier = trait, marker_name = markers[i], gxe = gxe_boolean, experiment_name = experiment)

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
	Take all the Markers and LOD scores above a certain value,
	
	The cutoff is given as input in the app
	'''

	cutoff_list = []
	# Iterate through the tuple_list and exclude all tuples with a
	#-logp value lower than the cutoff
	for info_tuple in tuple_list:#for each tuple (int, marker, lodscore)
		
		
		#Use abs() because there are negative -logp values present
		#The minus signs were added to discriminate between the parent alleles
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
	
	Also a single tuple is sent, however, it will be identified quickly
	as the maximum in the group.(Considering its the only one!)
	
	A list of 2 or more tuples will be fed to the function highest_tuple. 
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
		logp = [i[2] for i in region]
		max_logp = max(logp)
		logp_index = logp.index(max_logp)
		max_tuple = region[logp_index]
				
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
			#a physical position on a chromosome and is not related to
			#any marker and thus does not have a LOD score
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
	region_list : (trait, chromosome number, marker start position, marker end position)
	
	This information is used to collect lists of genes that are positioned
	on each region.(between the 2 physical start and end positions)
	[gene name, gene start, gene end]
	
	The tuple and gene lists will be placed together in a dictionary with
	the tuple as key and the list as value
	
	gene_dict = dict[(trait, chr nmbr, marker start pos, marker end pos)] = [gene name, gene start, gene end]
	
	NB: the chromosome number from the value is removed in func(normalize_object_data)
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
		normal_gene_list = normalize_object_data(remove_genes_behind)
		
		#add each region and list of genes to the dictionary gene_dict
		#with the tuple as a key and the gene list as a value
		gene_dict[region] = normal_gene_list
		
	return gene_dict
	
	
###############################################################################


def normalize_object_data(data_list):
	"""
	The extracted data from the database needs to be parsed 
	*data_list has the structure 
	[gene name,chromosome number, gene start, gene end]
	With data types that were derived from the model classes
	*The output should become 
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
