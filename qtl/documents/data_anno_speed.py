"""
This script is meant to test the speed of annotating a massive gene list

and possibly improve upon the speed


-edit: it worked:
Instead of reading a text file of 0.5 million lines
it reads a csv file with 20.000 rows

For a certain trait the annotation time decreased from 72 to 3 seconds!
"""

import time
import string
import csv



###################################################################	
###################################################################	


def read_once_text(filename):
	"""
	just a test...
	turns out it doesnt help
	"""
	
	tic = time.clock()
	
	fileobject = open(filename, 'r')
	data = fileobject.readlines()
	fileobject.close()
	
	how_many = len(data)
	
	
	toc = time.clock()
	stopwatch = toc-tic
	
	return data, how_many, stopwatch
	
	
def read_once_csv(filename):
	"""
	Read a csv file and place the contents in a dictionary
	"""
	
	tic = time.clock()
	
	#data = genfromtxt(filename, delimiter = ',')
	data_dict = {}
	
	with open(filename, 'rb') as csvobject:
		reader = csv.reader(csvobject)
		for row in reader:
			
			key = row[0]
			value = row[1:]
			
			data_dict[key] = value
		
	how_many = len(data_dict)
	
	toc = time.clock()
	stopwatch = toc-tic
	
	return data_dict, how_many, stopwatch


###################################################################	
###################################################################	
	

def annotate_from_text(data, gene_list):
	"""

	"""
	
	tic = time.clock()
	Gene_GO_list = []
	
	for gene_GO in data:

		s = string.split(gene_GO)
	
		#if the number of gene names in gene_list is more than zero
		if gene_list.count(s[0])>0:
		
			#only store the [genename, GOterm] part of the list
			info = s[:-1]
			Gene_GO_list.append(info)
			
		how_many = len(Gene_GO_list)
		
	toc = time.clock()
	stopwatch = toc-tic

	return how_many, stopwatch
	
	
	
def annotate_from_csv(data_dict, gene_list):
	"""
	This data comes from a csv file
	data is a dictionary
	each key is a gene name
	the corresponding value is a list of go terms
	"""
	tic = time.clock()
	
	gene_dict = {}	

	for gene in gene_list:
		
		if gene in data_dict.keys():
			
			#If the gene from the gene list is present in the data_dict
			#then add that gene as key to the new gene_dict with 
			#corresponding value
			gene_dict[gene] = data_dict[gene]
			


	how_many = len(gene_dict)

	toc = time.clock()
	stopwatch = toc-tic
	
	return gene_dict, how_many, stopwatch	
	

	
###################################################################	
###################################################################	


	
def print_file(trait, data):
	
	tic = time.clock()
	
	count_gene = 0
	
	for gene_GO in data:

		s = string.split(gene_GO)
		
		if trait in s:
			
			count_gene +=1
			
	toc = time.clock()
	stopwatch = toc-tic
			
	return count_gene, stopwatch


###################################################################	
###################################################################	


			
			
def populate_a_dict(data):
	"""
	Make a dictionary
	with genenames as key
	and a list of corresponding GO terms as value
	"""
	
	tic = time.clock()
	
	geneGO_dict = {}
	i = 1	
	
	for geneGO in data:
		
		s = string.split(geneGO)
		gene = s[0]
		GO = s[1]
		
		if gene in geneGO_dict.keys():
		#If the gene is already in the dictionary as a key than dont make 
		#a new key, all keys must be unique
		#Since the gene exists it must also have a value in the form of a
		#list containing GO annotations
		
			#Store the value in a temporary list
			temp_GO_list = geneGO_dict[gene]
			
		else:
			#If the gene is not yet in the disctionary as a key, than
			#use the gene name to make a new key.
			#The new key will need a value( list of GO annotations)
			#So make a new temporary empty list which will be populated
			#with its first GO term later on
			print i
			temp_GO_list = []
			i += 1
			
		
		#Add the GO term the new list
		temp_GO_list.append(GO)
		#Set the list of GO annotations as the value with the gene name
		#as key
		geneGO_dict[gene] = temp_GO_list
		
	toc = time.clock()
	stopwatch = toc-tic	
	
	
	#why not go ahead and write the data to a csv file immediately?
	#NO MODULARIZE!
	
	
	
	return geneGO_dict, stopwatch
	
	
	
def populate_a_list(gene_dict):
	"""
	Make a list
	With genenames on index position [0]
	And GO annotations on the positions after that
	"""
	tic = time.clock()
	
	geneGOarray = []
	for key, value in gene_dict.iteritems():
		
		
		geneGOlist = []
		geneGOlist.append(key)
		
		GO = gene_dict[key]
		#geneGOlist.append(GO)
		geneGOlist.extend(value)
			
		geneGOarray.append(geneGOlist)
		
	toc = time.clock()
	stopwatch = toc-tic	
	
	
	return geneGOarray, stopwatch	




#def populate_a_dict_from_csv(data):
	#"""
	
	#"""	
	
	#gene_list_array
	
	#for geneGO in data:
		
		#s = string.split(geneGO)
		#gene = s[0]
		#GO = s[1]




###################################################################	
###################################################################	

	

def write_file_for_gene(gene, golist):
	"""
	Take a gene (key) and a list of go annotations (value)
	And write them to a separate text file
	"""

	with open(gene, 'w') as fileobject:
		fileobject.write(gene)
		
		for go in golist:
			fileobject.write(go)
		

		
		
def write_dictionary_to_csv(gene_dict, filename):
	"""
	Write the dictionary to a csv file
	"""
	

	
	with open(filename, 'wb') as csvfile:
		writer = csv.DictWriter(csvfile, gene_dict)
		
		for key, value in gene_dict.iteritems():
			writer.writerow([key, value])
		
	
	
	
	
def write_list_to_csv(gene_array, filename):
	"""
	
	"""
	
	with open(filename, 'wb') as csvfile:
		
		for gene_list in gene_array:
		
			writer = csv.writer(csvfile)
			writer.writerow(gene_list)



def count_distinct_go(datadict):
	"""
	
	"""
	
	golist = []
	
	for key in datadict:
		
		for go in datadict[key]:
			
			if go not in golist:
				
				golist.append(go)
				
				
	return golist
	

###################################################################	
###################################################################	


	


	
if __name__ == "__main__":
	
	list10 = ['AT3G05327', 'AT3G05330', 'AT3G05340', 'AT3G05345', 'AT3G05350', 'AT3G05370', 'AT3G05380', 'AT3G05390', 'AT3G05400', 'AT3G05410']
	
	list100 = ['AT3G05327', 'AT3G05330', 'AT3G05340', 'AT3G05345', 'AT3G05350', 'AT3G05370', 'AT3G05380', 'AT3G05390', 'AT3G05400', 'AT3G05410', 'AT3G05420', 'AT3G05430', 'AT3G05440', 'AT3G05450', 'AT3G05460', 'AT3G05470', 'AT3G05480', 'AT3G05490', 'AT3G05500', 'AT3G05510', 'AT3G05520', 'AT3G05530', 'AT3G05540', 'AT3G05545', 'AT3G05550', 'AT3G05560', 'AT3G05570', 'AT3G05580', 'AT3G05590', 'AT3G05600', 'AT3G05602', 'AT3G05610', 'AT3G05620', 'AT3G05625', 'AT3G05630', 'AT3G05640', 'AT3G05650', 'AT3G05660', 'AT3G05670', 'AT3G05675', 'AT3G05680', 'AT3G05685', 'AT3G05690', 'AT3G05700', 'AT3G05710', 'AT3G05720', 'AT3G05725', 'AT3G05727', 'AT3G05730', 'AT3G05740', 'AT3G05741', 'AT3G05746', 'AT3G05750', 'AT3G05760', 'AT3G05770', 'AT3G05780', 'AT3G05790', 'AT3G05800', 'AT3G05810', 'AT3G05820', 'AT3G05830', 'AT3G05840', 'AT3G05858', 'AT3G05860', 'AT3G05880', 'AT3G05890', 'AT3G05900', 'AT3G05905', 'AT3G05910', 'AT3G05920', 'AT3G05930', 'AT3G05932', 'AT3G05935', 'AT3G05936', 'AT3G05937', 'AT3G05940', 'AT3G05950', 'AT3G05960', 'AT3G05970', 'AT3G05975', 'AT3G05980', 'AT3G05990', 'AT3G06000', 'AT3G06010', 'AT3G06020', 'AT3G06030', 'AT3G06035', 'AT3G06040', 'AT3G06050', 'AT3G06060', 'AT3G06070', 'AT3G06080', 'AT3G06090', 'AT3G06100', 'AT3G06110', 'AT3G06120', 'AT3G06125', 'AT3G06130', 'AT3G06140', 'AT3G06145']
	list1230 = ['AT3G05327', 'AT3G05330', 'AT3G05340', 'AT3G05345', 'AT3G05350', 'AT3G05370', 'AT3G05380', 'AT3G05390', 'AT3G05400', 'AT3G05410', 'AT3G05420', 'AT3G05430', 'AT3G05440', 'AT3G05450', 'AT3G05460', 'AT3G05470', 'AT3G05480', 'AT3G05490', 'AT3G05500', 'AT3G05510', 'AT3G05520', 'AT3G05530', 'AT3G05540', 'AT3G05545', 'AT3G05550', 'AT3G05560', 'AT3G05570', 'AT3G05580', 'AT3G05590', 'AT3G05600', 'AT3G05602', 'AT3G05610', 'AT3G05620', 'AT3G05625', 'AT3G05630', 'AT3G05640', 'AT3G05650', 'AT3G05660', 'AT3G05670', 'AT3G05675', 'AT3G05680', 'AT3G05685', 'AT3G05690', 'AT3G05700', 'AT3G05710', 'AT3G05720', 'AT3G05725', 'AT3G05727', 'AT3G05730', 'AT3G05740', 'AT3G05741', 'AT3G05746', 'AT3G05750', 'AT3G05760', 'AT3G05770', 'AT3G05780', 'AT3G05790', 'AT3G05800', 'AT3G05810', 'AT3G05820', 'AT3G05830', 'AT3G05840', 'AT3G05858', 'AT3G05860', 'AT3G05880', 'AT3G05890', 'AT3G05900', 'AT3G05905', 'AT3G05910', 'AT3G05920', 'AT3G05930', 'AT3G05932', 'AT3G05935', 'AT3G05936', 'AT3G05937', 'AT3G05940', 'AT3G05950', 'AT3G05960', 'AT3G05970', 'AT3G05975', 'AT3G05980', 'AT3G05990', 'AT3G06000', 'AT3G06010', 'AT3G06020', 'AT3G06030', 'AT3G06035', 'AT3G06040', 'AT3G06050', 'AT3G06060', 'AT3G06070', 'AT3G06080', 'AT3G06090', 'AT3G06100', 'AT3G06110', 'AT3G06120', 'AT3G06125', 'AT3G06130', 'AT3G06140', 'AT3G06145', 'AT3G06150', 'AT3G06160', 'AT3G06170', 'AT3G06180', 'AT3G06190', 'AT3G06200', 'AT3G06210', 'AT3G06220', 'AT3G06230', 'AT3G06240', 'AT3G06250', 'AT3G06260', 'AT3G06270', 'AT3G06280', 'AT3G06290', 'AT3G06300', 'AT3G06310', 'AT3G06320', 'AT3G06330', 'AT3G06340', 'AT3G06350', 'AT3G06360', 'AT3G06370', 'AT3G06380', 'AT3G06390', 'AT3G06400', 'AT3G06410', 'AT3G06420', 'AT3G06430', 'AT3G06435', 'AT3G06440', 'AT3G06450', 'AT3G06455', 'AT3G06460', 'AT3G06470', 'AT3G06480', 'AT3G06483', 'AT3G06490', 'AT3G06500', 'AT3G06510', 'AT3G06520', 'AT3G06530', 'AT3G06540', 'AT3G06545', 'AT3G06550', 'AT3G06560', 'AT3G06570', 'AT3G06580', 'AT3G06590', 'AT3G06600', 'AT3G06610', 'AT3G06620', 'AT3G06630', 'AT3G06640', 'AT3G06650', 'AT3G06660', 'AT3G06670', 'AT3G06680', 'AT3G06700', 'AT3G06710', 'AT3G06720', 'AT3G06730', 'AT3G06740', 'AT3G06750', 'AT3G06760', 'AT3G06770', 'AT3G06778', 'AT3G06780', 'AT3G06790', 'AT3G06810', 'AT3G06820', 'AT3G06830', 'AT3G06840', 'AT3G06850', 'AT3G06860', 'AT3G06868', 'AT3G06880', 'AT3G06890', 'AT3G06900', 'AT3G06910', 'AT3G06920', 'AT3G06930', 'AT3G06950', 'AT3G06960', 'AT3G06970', 'AT3G06980', 'AT3G06990', 'AT3G07000', 'AT3G07005', 'AT3G07010', 'AT3G07020', 'AT3G07025', 'AT3G07030', 'AT3G07040', 'AT3G07050', 'AT3G07060', 'AT3G07070', 'AT3G07080', 'AT3G07090', 'AT3G07100', 'AT3G07110', 'AT3G07120', 'AT3G07130', 'AT3G07140', 'AT3G07150', 'AT3G07160', 'AT3G07170', 'AT3G07180', 'AT3G07185', 'AT3G07190', 'AT3G07195', 'AT3G07200', 'AT3G07210', 'AT3G07215', 'AT3G07220', 'AT3G07230', 'AT3G07250', 'AT3G07260', 'AT3G07270', 'AT3G07273', 'AT3G07290', 'AT3G07300', 'AT3G07310', 'AT3G07320', 'AT3G07330', 'AT3G07340', 'AT3G07350', 'AT3G07360', 'AT3G07370', 'AT3G07390', 'AT3G07400', 'AT3G07410', 'AT3G07420', 'AT3G07425', 'AT3G07430', 'AT3G07440', 'AT3G07450', 'AT3G07460', 'AT3G07470', 'AT3G07480', 'AT3G07490', 'AT3G07500', 'AT3G07510', 'AT3G07520', 'AT3G07522', 'AT3G07525', 'AT3G07530', 'AT3G07540', 'AT3G07550', 'AT3G07560', 'AT3G07565', 'AT3G07568', 'AT3G07570', 'AT3G07580', 'AT3G07590', 'AT3G07600', 'AT3G07610', 'AT3G07620', 'AT3G07630', 'AT3G07640', 'AT3G07650', 'AT3G07660', 'AT3G07670', 'AT3G07680', 'AT3G07690', 'AT3G07700', 'AT3G07710', 'AT3G07720', 'AT3G07730', 'AT3G07740', 'AT3G07750', 'AT3G07760', 'AT3G07770', 'AT3G07780', 'AT3G07790', 'AT3G07800', 'AT3G07810', 'AT3G07820', 'AT3G07830', 'AT3G07840', 'AT3G07850', 'AT3G07860', 'AT3G07870', 'AT3G07880', 'AT3G07890', 'AT3G07900', 'AT3G07910', 'AT3G07920', 'AT3G07930', 'AT3G07940', 'AT3G07950', 'AT3G07960', 'AT3G07970', 'AT3G07980', 'AT3G07990', 'AT3G08000', 'AT3G08010', 'AT3G08020', 'AT3G08030', 'AT3G08040', 'AT3G08490', 'AT3G08500', 'AT3G08505', 'AT3G08510', 'AT3G08530', 'AT3G08550', 'AT3G08560', 'AT3G08570', 'AT3G08580', 'AT3G08590', 'AT3G08600', 'AT3G08610', 'AT3G08620', 'AT3G08630', 'AT3G08636', 'AT3G08640', 'AT3G08650', 'AT3G08660', 'AT3G08670', 'AT3G08680', 'AT3G08690', 'AT3G08700', 'AT3G08710', 'AT3G08720', 'AT3G08730', 'AT3G08740', 'AT3G08750', 'AT3G08760', 'AT3G08762', 'AT3G08770', 'AT3G08780', 'AT3G08790', 'AT3G08800', 'AT3G08810', 'AT3G08820', 'AT3G08840', 'AT3G08850', 'AT3G08860', 'AT3G08870', 'AT3G08880', 'AT3G08890', 'AT3G08900', 'AT3G08910', 'AT3G08920', 'AT3G08930', 'AT3G08940', 'AT3G08943', 'AT3G08947', 'AT3G08950', 'AT3G08960', 'AT3G08970', 'AT3G08980', 'AT3G08990', 'AT3G09000', 'AT3G09010', 'AT3G09020', 'AT3G09030', 'AT3G09032', 'AT3G09035', 'AT3G09040', 'AT3G09050', 'AT3G09060', 'AT3G09070', 'AT3G09080', 'AT3G09085', 'AT3G09090', 'AT3G09100', 'AT3G09110', 'AT3G09120', 'AT3G09130', 'AT3G09140', 'AT3G09150', 'AT3G09160', 'AT3G09162', 'AT3G09180', 'AT3G09190', 'AT3G09200', 'AT3G09210', 'AT3G09220', 'AT3G09230', 'AT3G09240', 'AT3G09250', 'AT3G09260', 'AT3G09270', 'AT3G09280', 'AT3G09290', 'AT3G09300', 'AT3G09310', 'AT3G09320', 'AT3G09330', 'AT3G09340', 'AT3G09350', 'AT3G09360', 'AT3G09370', 'AT3G09380', 'AT3G09390', 'AT3G09400', 'AT3G09405', 'AT3G09410', 'AT3G09430', 'AT3G09440', 'AT3G09450', 'AT3G09470', 'AT3G09480', 'AT3G09490', 'AT3G09500', 'AT3G09520', 'AT3G09530', 'AT3G09540', 'AT3G09550', 'AT3G09560', 'AT3G09570', 'AT3G09580', 'AT3G09590', 'AT3G09600', 'AT3G09620', 'AT3G09630', 'AT3G09640', 'AT3G09650', 'AT3G09660', 'AT3G09670', 'AT3G09680', 'AT3G09690', 'AT3G09700', 'AT3G09710', 'AT3G09720', 'AT3G09730', 'AT3G09735', 'AT3G09740', 'AT3G09750', 'AT3G09760', 'AT3G09770', 'AT3G09780', 'AT3G09790', 'AT3G09800', 'AT3G09810', 'AT3G09820', 'AT3G09830', 'AT3G09840', 'AT3G09850', 'AT3G09860', 'AT3G09863', 'AT3G09870', 'AT3G09880', 'AT3G09890', 'AT3G09900', 'AT3G09910', 'AT3G09920', 'AT3G09922', 'AT3G09925', 'AT3G09930', 'AT3G09940', 'AT3G09950', 'AT3G09960', 'AT3G09970', 'AT3G09975', 'AT3G09980', 'AT3G09990', 'AT3G10000', 'AT3G10010', 'AT3G10015', 'AT3G10020', 'AT3G10030', 'AT3G10040', 'AT3G10050', 'AT3G10060', 'AT3G10070', 'AT3G10080', 'AT3G10116', 'AT3G10120', 'AT3G10130', 'AT3G10140', 'AT3G10150', 'AT3G10160', 'AT3G10180', 'AT3G10190', 'AT3G10195', 'AT3G10200', 'AT3G10210', 'AT3G10220', 'AT3G10230', 'AT3G10240', 'AT3G10250', 'AT3G10260', 'AT3G10270', 'AT3G10280', 'AT3G10290', 'AT3G10300', 'AT3G10310', 'AT3G10320', 'AT3G10330', 'AT3G10340', 'AT3G10350', 'AT3G10360', 'AT3G10370', 'AT3G10380', 'AT3G10390', 'AT3G10400', 'AT3G10405', 'AT3G10410', 'AT3G10420', 'AT3G10430', 'AT3G10439', 'AT3G10440', 'AT3G10450', 'AT3G10460', 'AT3G10470', 'AT3G10480', 'AT3G10490', 'AT3G10500', 'AT3G10510', 'AT3G10520', 'AT3G10525', 'AT3G10526', 'AT3G10530', 'AT3G10540', 'AT3G10550', 'AT3G10560', 'AT3G10570', 'AT3G10572', 'AT3G10580', 'AT3G10585', 'AT3G10590', 'AT3G10595', 'AT3G10600', 'AT3G10610', 'AT3G10620', 'AT3G10630', 'AT3G10640', 'AT3G10650', 'AT3G10660', 'AT3G10670', 'AT3G10680', 'AT3G10690', 'AT3G10700', 'AT3G10710', 'AT3G10720', 'AT3G10730', 'AT3G10740', 'AT3G10745', 'AT3G10750', 'AT3G10760', 'AT3G10770', 'AT3G10780', 'AT3G10790', 'AT3G10800', 'AT3G10810', 'AT3G10815', 'AT3G10820', 'AT3G10830', 'AT3G10840', 'AT3G10845', 'AT3G10850', 'AT3G10860', 'AT3G10870', 'AT3G10880', 'AT3G10890', 'AT3G10900', 'AT3G10910', 'AT3G10915', 'AT3G10920', 'AT3G10930', 'AT3G10940', 'AT3G10950', 'AT3G10960', 'AT3G10970', 'AT3G10974', 'AT3G10980', 'AT3G10985', 'AT3G10986', 'AT3G10990', 'AT3G11000', 'AT3G11010', 'AT3G11020', 'AT3G11030', 'AT3G11040', 'AT3G11050', 'AT3G11070', 'AT3G11080', 'AT3G11090', 'AT3G11100', 'AT3G11110', 'AT3G11120', 'AT3G11130', 'AT3G11150', 'AT3G11160', 'AT3G11165', 'AT3G11170', 'AT3G11180', 'AT3G11200', 'AT3G11210', 'AT3G11220', 'AT3G11230', 'AT3G11240', 'AT3G11250', 'AT3G11260', 'AT3G11270', 'AT3G11280', 'AT3G11290', 'AT3G11300', 'AT3G11310', 'AT3G11320', 'AT3G11325', 'AT3G11330', 'AT3G11340', 'AT3G11350', 'AT3G11370', 'AT3G11380', 'AT3G11385', 'AT3G11390', 'AT3G11397', 'AT3G11400', 'AT3G11402', 'AT3G11405', 'AT3G11410', 'AT3G11420', 'AT3G11430', 'AT3G11435', 'AT3G11440', 'AT3G11450', 'AT3G11460', 'AT3G11470', 'AT3G11480', 'AT3G11490', 'AT3G11500', 'AT3G11510', 'AT3G11520', 'AT3G11530', 'AT3G11540', 'AT3G11550', 'AT3G11560', 'AT3G11570', 'AT3G11580', 'AT3G11590', 'AT3G11591', 'AT3G11600', 'AT3G11620', 'AT3G11630', 'AT3G11640', 'AT3G11650', 'AT3G11660', 'AT3G11670', 'AT3G11680', 'AT3G11690', 'AT3G11700', 'AT3G11710', 'AT3G11720', 'AT3G11730', 'AT3G11740', 'AT3G11745', 'AT3G11750', 'AT3G11760', 'AT3G11770', 'AT3G11773', 'AT3G11780', 'AT3G11800', 'AT3G11810', 'AT3G11820', 'AT3G11825', 'AT3G11830', 'AT3G11840', 'AT3G11850', 'AT3G11860', 'AT3G11870', 'AT3G11880', 'AT3G11890', 'AT3G11900', 'AT3G11910', 'AT3G11920', 'AT3G11930', 'AT3G11940', 'AT3G11945', 'AT3G11950', 'AT3G11960', 'AT3G11964', 'AT3G11980', 'AT3G12000', 'AT3G12010', 'AT3G12020', 'AT3G12030', 'AT3G12040', 'AT3G12050', 'AT3G12060', 'AT3G12070', 'AT3G12080', 'AT3G12090', 'AT3G12100', 'AT3G12110', 'AT3G12120', 'AT3G12130', 'AT3G12140', 'AT3G12145', 'AT3G12150', 'AT3G12160', 'AT3G12170', 'AT3G12180', 'AT3G12190', 'AT3G12200', 'AT3G12203', 'AT3G12210', 'AT3G12220', 'AT3G12230', 'AT3G12240', 'AT3G12250', 'AT3G12260', 'AT3G12270', 'AT3G12280', 'AT3G12290', 'AT3G12300', 'AT3G12320', 'AT3G12340', 'AT3G12345', 'AT3G12350', 'AT3G12360', 'AT3G12370', 'AT3G12380', 'AT3G12390', 'AT3G12400', 'AT3G12410', 'AT3G12420', 'AT3G12430', 'AT3G12440', 'AT3G12460', 'AT3G12470', 'AT3G12480', 'AT3G12490', 'AT3G12500', 'AT3G12510', 'AT3G12520', 'AT3G12530', 'AT3G12540', 'AT3G12545', 'AT3G12550', 'AT3G12560', 'AT3G12570', 'AT3G12580', 'AT3G12587', 'AT3G12590', 'AT3G12600', 'AT3G12610', 'AT3G12620', 'AT3G12630', 'AT3G12640', 'AT3G12650', 'AT3G12660', 'AT3G12670', 'AT3G12680', 'AT3G12685', 'AT3G12690', 'AT3G12700', 'AT3G12710', 'AT3G12720', 'AT3G12730', 'AT3G12740', 'AT3G12750', 'AT3G12760', 'AT3G12770', 'AT3G12775', 'AT3G12780', 'AT3G12800', 'AT3G12810', 'AT3G12820', 'AT3G12830', 'AT3G12840', 'AT3G12850', 'AT3G12860', 'AT3G12870', 'AT3G12880', 'AT3G12890', 'AT3G12900', 'AT3G12910', 'AT3G12915', 'AT3G12920', 'AT3G12930', 'AT3G12940', 'AT3G12950', 'AT3G12955', 'AT3G12960', 'AT3G12970', 'AT3G12977', 'AT3G12980', 'AT3G12990', 'AT3G13000', 'AT3G13010', 'AT3G13020', 'AT3G13030', 'AT3G13040', 'AT3G13050', 'AT3G13060', 'AT3G13062', 'AT3G13065', 'AT3G13070', 'AT3G13080', 'AT3G13090', 'AT3G13100', 'AT3G13110', 'AT3G13120', 'AT3G13130', 'AT3G13140', 'AT3G13150', 'AT3G13160', 'AT3G13170', 'AT3G13175', 'AT3G13180', 'AT3G13190', 'AT3G13200', 'AT3G13210', 'AT3G13220', 'AT3G13222', 'AT3G13224', 'AT3G13225', 'AT3G13226', 'AT3G13227', 'AT3G13228', 'AT3G13229', 'AT3G13230', 'AT3G13235', 'AT3G13240', 'AT3G13275', 'AT3G13277', 'AT3G13280', 'AT3G13290', 'AT3G13300', 'AT3G13310', 'AT3G13320', 'AT3G13330', 'AT3G13340', 'AT3G13350', 'AT3G13360', 'AT3G13370', 'AT3G13380', 'AT3G13390', 'AT3G13400', 'AT3G13404', 'AT3G13405', 'AT3G13410', 'AT3G13420', 'AT3G13430', 'AT3G13432', 'AT3G13433', 'AT3G13435', 'AT3G13437', 'AT3G13440', 'AT3G13445', 'AT3G13450', 'AT3G13460', 'AT3G13470', 'AT3G13480', 'AT3G13490', 'AT3G13500', 'AT3G13510', 'AT3G13520', 'AT3G13525', 'AT3G13530', 'AT3G13540', 'AT3G13550', 'AT3G13560', 'AT3G13570', 'AT3G13580', 'AT3G13590', 'AT3G13594', 'AT3G13600', 'AT3G13610', 'AT3G13620', 'AT3G13630', 'AT3G13640', 'AT3G13650', 'AT3G13660', 'AT3G13662', 'AT3G13670', 'AT3G13672', 'AT3G13674', 'AT3G13677', 'AT3G13680', 'AT3G13682', 'AT3G13690', 'AT3G13710', 'AT3G13720', 'AT3G13724', 'AT3G13730', 'AT3G13740', 'AT3G13750', 'AT3G13760', 'AT3G13770', 'AT3G13772', 'AT3G13780', 'AT3G13782', 'AT3G13784', 'AT3G13790', 'AT3G13800', 'AT3G13810', 'AT3G13820', 'AT3G13830', 'AT3G13840', 'AT3G13845', 'AT3G13850', 'AT3G13855', 'AT3G13860', 'AT3G13870', 'AT3G13880', 'AT3G13882', 'AT3G13890', 'AT3G13900', 'AT3G13910', 'AT3G13920', 'AT3G13930', 'AT3G13940', 'AT3G13950', 'AT3G13960', 'AT3G13970', 'AT3G13980', 'AT3G13990', 'AT3G14000', 'AT3G14010', 'AT3G14020', 'AT3G14030', 'AT3G14040', 'AT3G14050', 'AT3G14060', 'AT3G14067', 'AT3G14070', 'AT3G14075', 'AT3G14080', 'AT3G14090', 'AT3G14100', 'AT3G14110', 'AT3G14120', 'AT3G14130', 'AT3G14140', 'AT3G14150', 'AT3G14160', 'AT3G14170', 'AT3G14172', 'AT3G14180', 'AT3G14185', 'AT3G14190', 'AT3G14200', 'AT3G14205', 'AT3G14210', 'AT3G14220', 'AT3G14225', 'AT3G14230', 'AT3G14240', 'AT3G14250', 'AT3G14260', 'AT3G14270', 'AT3G14280', 'AT3G14290', 'AT3G14300', 'AT3G14310', 'AT3G14320', 'AT3G14330', 'AT3G14350', 'AT3G14360', 'AT3G14362', 'AT3G14370', 'AT3G14380', 'AT3G14385', 'AT3G14390', 'AT3G14395', 'AT3G14400', 'AT3G14410', 'AT3G14415', 'AT3G14420', 'AT3G14430', 'AT3G14431', 'AT3G14440', 'AT3G14450', 'AT3G14452', 'AT3G14460', 'AT3G14470', 'AT3G14480', 'AT3G14490', 'AT3G14510', 'AT3G14520', 'AT3G14530', 'AT3G14540', 'AT3G14550', 'AT3G14560', 'AT3G14570', 'AT3G14580', 'AT3G14590', 'AT3G14595', 'AT3G14600', 'AT3G14610', 'AT3G14620', 'AT3G14630', 'AT3G14640', 'AT3G14650', 'AT3G14660', 'AT3G14670', 'AT3G14680', 'AT3G14690', 'AT3G14700', 'AT3G14710', 'AT3G14720', 'AT3G14730', 'AT3G14735', 'AT3G14740', 'AT3G14750', 'AT3G14760', 'AT3G14770', 'AT3G14780', 'AT3G14790', 'AT3G14810', 'AT3G14820', 'AT3G14830', 'AT3G14840', 'AT3G14850', 'AT3G14860', 'AT3G14870', 'AT3G14880', 'AT3G14890', 'AT3G14900', 'AT3G14910', 'AT3G14920', 'AT3G14930', 'AT3G14940', 'AT3G14950', 'AT3G14960', 'AT3G14970', 'AT3G14980', 'AT3G14990', 'AT3G15000', 'AT3G15010', 'AT3G15020', 'AT3G15030', 'AT3G15040', 'AT3G15050', 'AT3G15060', 'AT3G15070', 'AT3G15080', 'AT3G15090', 'AT3G15095', 'AT3G15110', 'AT3G15111', 'AT3G15115', 'AT3G15120', 'AT3G15130', 'AT3G15140', 'AT3G15150', 'AT3G15160', 'AT3G15170', 'AT3G15180', 'AT3G15190', 'AT3G15200', 'AT3G15210', 'AT3G15220', 'AT3G15240', 'AT3G15250', 'AT3G15251', 'AT3G15260', 'AT3G15280', 'AT3G15290', 'AT3G15300', 'AT3G15340', 'AT3G15350', 'AT3G15351', 'AT3G15352', 'AT3G15353', 'AT3G15354', 'AT3G15355', 'AT3G15356', 'AT3G15357', 'AT3G15358', 'AT3G15359', 'AT3G15360', 'AT3G15370', 'AT3G15380', 'AT3G15390', 'AT3G15395', 'AT3G15400', 'AT3G15410', 'AT3G15420', 'AT3G15430', 'AT3G15440', 'AT3G15450', 'AT3G15460', 'AT3G15470', 'AT3G15480', 'AT3G15490', 'AT3G15500', 'AT3G15510', 'AT3G15518', 'AT3G15520', 'AT3G15530', 'AT3G15534', 'AT3G15540', 'AT3G15548', 'AT3G15550', 'AT3G15570', 'AT3G15578', 'AT3G15580', 'AT3G15590', 'AT3G15605', 'AT3G15610', 'AT3G15620', 'AT3G15630', 'AT3G15640', 'AT3G15650', 'AT3G15660', 'AT3G15670', 'AT3G15680', 'AT3G15690', 'AT3G15700', 'AT3G15710', 'AT3G15720', 'AT3G15730', 'AT3G15740', 'AT3G15760', 'AT3G15770', 'AT3G15780', 'AT3G15790', 'AT3G15800', 'AT3G15810', 'AT3G15820', 'AT3G15830', 'AT3G15840', 'AT3G15850', 'AT3G15860', 'AT3G15870', 'AT3G15880', 'AT3G15890', 'AT3G15900', 'AT3G15910', 'AT3G15920', 'AT3G15930', 'AT3G15940', 'AT3G15950', 'AT3G15960', 'AT3G15970', 'AT3G15980', 'AT3G15990', 'AT3G16000', 'AT3G16010', 'AT3G16020', 'AT3G16030', 'AT3G16040', 'AT3G16050', 'AT3G16060', 'AT3G16070', 'AT3G16080', 'AT3G16090', 'AT3G16100', 'AT3G16110', 'AT3G16117', 'AT3G16120', 'AT3G16130', 'AT3G16140', 'AT3G16150', 'AT3G16160', 'AT3G16170', 'AT3G16175', 'AT3G16180', 'AT3G16190', 'AT3G16200', 'AT3G16210', 'AT3G16220', 'AT3G16230', 'AT3G16240', 'AT3G16250', 'AT3G16260', 'AT3G16270', 'AT3G16280', 'AT3G16290', 'AT3G16300', 'AT3G16310', 'AT3G16320', 'AT3G16330', 'AT3G16340', 'AT3G16350', 'AT3G16360', 'AT3G16370', 'AT3G16380', 'AT3G16390', 'AT3G16400', 'AT3G16410', 'AT3G16420', 'AT3G16430', 'AT3G16440', 'AT3G16450', 'AT3G16460', 'AT3G16470', 'AT3G16480', 'AT3G16490', 'AT3G16500', 'AT3G16510', 'AT3G16520', 'AT3G16530', 'AT3G16540', 'AT3G16550', 'AT3G16555', 'AT3G16560', 'AT3G16565', 'AT3G16570', 'AT3G16580', 'AT3G16590', 'AT3G16600', 'AT3G16610', 'AT3G16620', 'AT3G16630', 'AT3G16640', 'AT3G16650', 'AT3G16660', 'AT3G16670', 'AT3G16680', 'AT3G16690', 'AT3G16700', 'AT3G16710', 'AT3G16712', 'AT3G16720', 'AT3G16730', 'AT3G16740', 'AT3G16750', 'AT3G16760', 'AT3G16770', 'AT3G16780', 'AT3G16785', 'AT3G16800', 'AT3G16810', 'AT3G16820', 'AT3G16830', 'AT3G16840', 'AT3G16850', 'AT3G16851', 'AT3G36659', 'AT3G66652', 'AT3G66654', 'AT3G66656', 'AT3G66658']
	
	#dicttest = {"aaa":list10 , "bbb":list100, "ccc":list1230 ,"ddd":list100 ,"eee":list10}
	
	#filename = 'ArabiGOannotations.out'
	
	#csv_filename1 = "tester1.csv"
	csv_filename2 = 'At_geneGOlists.csv'

	
	#data, amount, timed = read_once(filename)
	#print "The data contains %d lines" % amount
	#print "Time to fill data is %f seconds" % timed
	#print "################################################"
	
	
	#gene_dict, timed_dict = populate_a_dict(data)
	#print "The dictionary of genes contains %d keys" % len(gene_dict)
	#print "Time to fill the dictionary is %f seconds" % timed_dict
	#print "################################################"
	
	
	#gene_array, timed_list = populate_a_list(gene_dict)
	#print "The array of genes contains %d lists" % len(gene_array)
	#print "Time to fill the array is %f seconds" % timed_list
	#print "################################################"
	
	
	#write_list_to_csv(gene_array, csv_filename2)
	
	
	#write_dictionary_to_csv(dicttest, csv_filename1)
	

	#write_list_to_csv(garray, csv_filename2)
	
	
######################################################
######################################################
######################################################
######################################################

	datadct, timed, stopwatch = read_once_csv(csv_filename2)
	
	
	golist = count_distinct_go(datadct)
	
	print "number of distinct go terms is %d" % len(golist)
	#= 626
	
	
	#t = "AT1G01060"
	#n, timed_n = print_file(t, data)
	
	
	#print "%s is mentioned %d times in this file" % (t, n)
	#print "Time to read the file is %f seconds" % timed
	#print "################################################"	
	
	
	#print "There are %d lines in this file" % amount
	#print "Time to read the file is %f seconds" % timed
	#print "################################################"
	
	#gene_dict10, amount_10, timed_10 = annotate_from_csv(datadct, list10)
	
	#print 
	#print "Number of genes in the gene list is %d" % len(list10)
	#print "There are %d annotations in this file" % amount_10
	#print "Time to annotate the genes is %f seconds" % timed_10
	#print "################################################"
	
	#gene_dict100, amount_100, timed_100 = annotate_from_csv(datadct, list100)
	
	#print "Number of genes in the gene list is %d" % len(list100)
	#print "There are %d annotations in this file" % amount_100
	#print "Time to annotate the genes is %f seconds" % timed_100
	#print "################################################"
	
	#gene_dict1230, amount_1230, timed_1230 = annotate_from_csv(datadct, list1230)
	
	#print "Number of genes in the gene list is %d" % len(list1230)
	#print "There are %d annotations in this file" % amount_1230
	#print "Time to annotate the genes is %f seconds" % timed_1230
	#print "################################################"



	#print "gene_dict10 = ", gene_dict10






