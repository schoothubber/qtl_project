import sys
import os
import csv
from collections import Counter
from sortedcontainers import SortedDict
from scipy import stats
from numpy import array, empty 

#from qtl.models import Gene,Marker,LOD,Experiment,ExperimentMarker


sys.path.append('/home/wouter/xenv_dj1.6/QTL')
#sys.path.append('/mnt/geninf15/prog/www/django/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')


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

	data_dict : dict[gene] = [go_list]
	gene_list : [list of genes from the eQTLs]
	
	The output are two dictionaries for genes inside and outside the eQTL's
	
	gene_inside_qtl_dict : dict[gene] = [go_list]
	gene_outside_qtl_dict : dict[gene] = [go_list]
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
	


def read_go_info_once_csv(filename):
	"""
	Read a csv file and place the contents in a dictionary
	
	The csv file contains rows of data
	In each row the first index holds a gene name
	The rest of the row holds the GO annotations
	
	data_dict : dict[go term] = [name, definition]
	"""
	
	
	data_dict = {}
	
	with open(filename, 'rb') as csvobject:
		reader = csv.reader(csvobject)
		for row in reader:
			
			#key is the go term
			key = row[0]
			#value is a name and a definition
			value = [row[1], row[3]]
			
			data_dict[key] = value
	
	
	return data_dict



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




def flatten_array(gene_in_qtl_dict):
	"""
	Make go_array from dict values
	The values are lists of go terms
	
	Note that individual GO terms can be present in the list multiple 
	times. This is exactly the point
	"""
	go_array = []
	
	for key, value in gene_in_qtl_dict.iteritems():
		
		go_array += value
	
	return go_array



	
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
	
	The output is a go term and the corresponding RPV
	
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
		
		#[i, go, fu_pvalue]
		
	
	return fisher_output
	
	
def extract_significant_result_fish(datainput, alpha):
	"""
	Only allow the results with a pvalue equal or below alpha to be returned
	Because those are to be considered significant results
	"""
	
	significant_output = []
	
	#info = [i, go, fu_pvalue]
	
	
	for info in datainput:
		if info[2] <= alpha:
			
		
			significant_output.append(info)
			
	
	return significant_output
	
def extract_significant_result_mult(datainput, alpha):
	"""
	Only allow the results with a pvalue equal or below alpha to be returned
	Because those are to be considered significant results
	"""
	
	significant_output = []
	
	#info = [i, go, fu_pvalue]
	
	
	for info in datainput:
		if info[3] <= alpha:
			
		
			significant_output.append(info)
			
	
	return significant_output
	

def correct_pvalues_for_multiple_testing(data):                
	"""                                                                                                   
	consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
	"""

	pvalues = [info[2] for info in data]

	pvalues = array(pvalues) 
	n = float(pvalues.shape[0])                                                                           
	new_pvalues = empty(n)
                                                                                                                       
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ] 
	                                     
	values.sort()
	values.reverse()                                                                                  
	new_values = []
	
	for i, vals in enumerate(values):                                                                 
		rank = n - i
		pvalue, index = vals                                                                          
		new_values.append((n/rank) * pvalue)
		                                                          
	for i in xrange(0, int(n)-1):  
		if new_values[i] < new_values[i+1]:                                                           
			new_values[i+1] = new_values[i] 
			                                                          
	for i, vals in enumerate(values):
		pvalue, index = vals
		new_pvalues[index] = new_values[i] 
	
	adjusted_output = []
	for i in range(0, len(data)):
		
		q_value = round(new_pvalues[i], 7)
		individual_adjusted = [data[i][0], data[i][1], data[i][2] , q_value]
		
		adjusted_output.append(individual_adjusted)
		
		                                                                                                                 
	return adjusted_output



	
###############################################################################
############################Post Processing#####################################
####################AaltJan has some demands!###################################
###############################################################################
	




def make_qtl_go_dict(gene_in_qtl_dict, gene_dict):
	"""
	Make a new dictionary with QTLs as keys.
	replace the values of genelists with corresponding golists
	Resulting in a dictionary of QTLs with a list of the found go terms
	for each QTL
	
	gene_in_qtl_dict : dict[gene] = [go_list]
	gene_dict : dict[(trait, chr nmbr, marker start pos, marker end pos)] = [gene name, gene start, gene end]
	
	qtl_go_dict : dict[(trait, chr nmbr, marker start pos, marker end pos)] = [go_list] 
	
	NB: each go list from values in qtl_go_dict do not contain unique go terms!!
	"""
	
	qtl_go_dict = {}
	
	#Each of these qtls contain a gene_list which will be iterated through
	for qtl in gene_dict:
		
		gene_loc_list = gene_dict[qtl]
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
	
	enriched_golist : [i, go, fu_p_value[0], adjusted(fu_p_value[0]) ]
	qtl_go_dict : dict[(trait, chr nmbr, marker start pos, marker end pos)] = [go_list]
	
	approved_golist_A : [i, go, fu_p_value[0], adjusted(fu_p_value[0]), go_fraction_score]
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
	For each GO term count the number of genes that have this GO term
	for the entire genome!
	
	Divide this number through the total number of genes
	(Question: total number of known genes or total number of genes used)
	
	approved_golist_A : [i, go, fu_p_value[0], adjusted(fu_p_value[0]), go_fraction_score]
	gene_in_qtl_dict : dict[gene] = [go_list]
	gene_out_qtl_dict : dict[gene] = [go_list]
	
	
	go_gene_dict_full : dict[i, go, fu_p_value, adj(fu_p_val, go_frac_scA, go_frac_scB]) = [gene list]

	
	"""
	
	#Count for how many genes in the QTL a go term from golist is predicted
	#Need approved_golist_A + gene_in_qtl_dict
	
	approved_golist_B = []
	n_genes_in_qtl = len(gene_in_qtl_dict)
	n_genes_out_qtl = len(gene_out_qtl_dict)	
	total_number_genes = float(n_genes_in_qtl + n_genes_out_qtl)
	#print "total genes:", total_number_genes
	
	#Inherit the SortedDict which will enable the keys to remain sorted
	go_gene_dict_full = SortedDict()
	
	for go_info in approved_golist_A:
		
		go = go_info[1]
		go_count = 0.0
		gene_list_in = []

		
		for gene in gene_in_qtl_dict:
			
			if go in gene_in_qtl_dict[gene]:
				go_count += 1.0
				#only add a gene if its found in a eQTL
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

			
	return go_gene_dict_full
