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
	
