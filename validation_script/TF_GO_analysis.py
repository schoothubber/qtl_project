import sys
import os

from sortedcontainers import SortedDict

from qtl.compare_with_literature import parse_data

#The following lines enable this script to contact django
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')


def read_data(fn):
	"""
	Eat and spit data from a file
	"""
	with open(fn, 'r') as fo:
		data = fo.readlines()
	return data
	
	
def write_location_data(fn, result):
	"""
	Document the results in a text file
	result is a dictionary
	"""
	with open(fn, 'w') as fo:
		for key, value in result.iteritems():
			fo.write('trait: %s'%key[0])
			fo.write('\n')
			fo.write('location: %s'%key[1])
			fo.write('\n')
			fo.write('TF: ')
			for tf in value:
				fo.write('%s '%tf)
			fo.write('\n\n')
			
def write_data(fn, result):
	"""
	Document the results in a text file
	result is a dictionary
	"""
	with open(fn, 'w') as fo:
		for key, value in result.iteritems():
			fo.write('trait: %s'%key)
			fo.write('\n')
			fo.write('TF: ')
			for tf in value:
				fo.write('%s '%tf)
			fo.write('\n\n')
			
def write_data2(fn, result):
	"""
	Document the results in a text file
	result is an array
	"""
	with open(fn, 'w') as fo:
		for cell in result:
			for item in cell:
				fo.write("%s "%item)
			fo.write('\n\n')
		

def get_go_terms(fn):
	"""
	Retrieve all connections from go-basic.obo
	"""
	
	data = read_data(fn)
	
	is_adict = SortedDict()
	is_alist = []
	
	for line in data:
		if line.startswith("id: GO:"):
			go = line[4:14]
		if line.startswith("is_a: GO:"):
			is_ago = line[6:16]
			
			if not go in is_adict:
				is_alist = []
				is_alist.append(is_ago)
				is_adict[go] = is_alist
			else:
				is_atemp = is_adict[go]
				is_atemp.append(is_ago)
				is_adict[go] = is_atemp
				
	#print is_adict
				
	return is_adict
				

def check_TF_results_cistrans(fn):
	"""
	Divide the found regulators into cis and trans
	This happens per trait
	"""
	
	data = read_data(fn)
	
	tf_dict = {}
	
	for line in data:
		
		#Read and Arrange the data:
		if line.startswith("TGR:"):
			info = line[5:].split()
			
			trait = info[0]
			go = info[1]
			region = info[2]
				
		if line.startswith("AT"):
			tf = line.split()
			tf_set = set(tf)
			
			#Start to manipulate:
			trait_chrom = trait[2]
			eQTL_chrom = region[0]
			
			#Divide the regulators in cis and trans for each trait:
			#-------------------------------------------------------
			#NB:Take care in the selection of cis eQTLs.
			#The eQTL intervals contain alot of genes
			#Therefor there could be overlap with local-trans eQTLs
			#-------------------------------------------------------
			if trait_chrom == eQTL_chrom:#it shall be cis
				cis_tup = (trait, "cis")#the key for cis eQTLs
				
				if cis_tup not in tf_dict:
					tf_dict[cis_tup] = tf_set
				else:
					cis_temp = tf_dict[cis_tup]
					cistf_union = cis_temp | tf_set #a union of sets
					tf_dict[cis_tup] = cistf_union
				
			else:#it shall be trans
				trans_tup = (trait, "trans")#the Key for Trans eQTLs
				
				if trans_tup not in tf_dict:
					tf_dict[trans_tup] = tf_set
				else:
					trans_temp = tf_dict[trans_tup]
					transtf_union = trans_temp | tf_set #a union of sets
					tf_dict[trans_tup] = transtf_union			

	#print tf_dict
	return tf_dict


def make_TF_dict(fn):
	"""
	Collect all found regulators per trait
	For comparison with confirmed data
	"""
	data = read_data(fn)
	
	tf_dict = {}
	
	for line in data:
		
		#Read and Arrange the data:
		if line.startswith("TGR:"):
			info = line[5:].split()
			
			trait = info[0]
			go = info[1]
			region = info[2]
				
		if line.startswith("AT"):
			tf = line.split()
			tf_set = set(tf)
		
			#Collect the data:
			if trait not in tf_dict:
				tf_dict[trait] = tf_set
			else:
				temp_reg = tf_dict[trait]
				tf_union = temp_reg | tf_set #a union of sets
				tf_dict[trait] = tf_union
			
	return tf_dict
			

def sort_dict_and_sets(adict):
	"""
	Take an unsorted dictionary and give back a sorted dictionary
	The keys will be sorted and the values will be sorted
	"""
	
	bdict = SortedDict()
	for key, value in adict.iteritems():
		bdict[key] = sorted(value)
	return bdict
		
	


		
"""
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/TF_GO_analysis.py')
"""



#filename = "qtl/documents/go-basic.obo"
#get_go_terms(filename)


exp_list = ['Ligterink_2014','Keurentjes_2007','Snoek_2012']


dataset = exp_list[0]
folder = "enrichment_results/"
filein1 = "%sTF_results_%s.txt"%(folder, dataset)
filein2 = 'raw_data/AtRegNet.txt'
fileout1 = "%sTF_per_trait_%s.txt"%(folder, dataset)
fileout2 = "%sTF_validation_%s.txt"%(folder, dataset)




#Write the results to a file:
#print "Saving file..."
#write_data(fileout, sorted_tfdict)





tfdict = make_TF_dict(filein1)
#sorted_tfdict = sort_dict_and_sets(tfdict)
	

#note:
	#overlapping genes
	#missing genes
	#number of overlapping genes
	#number of actual genes
	
confirm_data = []
overlap_dict = {}
missing_dict = {}
for gene, tf in confirmed_TF_dict.iteritems():
	if gene in tfdict:
		#----------------------------
		#Make sets
		tf_check = tfdict[gene]
		tf_check_set = set(tf_check)
		tf_set = set(tf)
		#----------------------------
		#The actual genes:
		nr_confirmed_tf = len(tf_set)
		#The overlapping genes:
		tf_overlap = tf_check_set & tf_set
		nr_tf_overlap = len(tf_overlap)
		#The missing genes:
		tf_missing = tf_check_set - tf_set
		#Percentage confirmed
		if nr_tf_overlap != 0:
			perc = float(nr_tf_overlap/nr_confirmed_tf)
		else:
			perc = 0
		
		confirm_cell = [gene, nr_tf_overlap, nr_confirmed_tf, perc] 
		confirm_data.append(confirm_cell)
		
		overlap_dict[gene] = tf_overlap
		missing_dict[gene] = tf_missing
		
#for k,v in overlap_dict.iteritems():
	#print "k: %s"%k
	##for val in v:
	#print "v: %s"%v
		
		
#write_data2(fileout2, confirm_data)			
	
#t = 'AT5G44030'

#if t in confirmed_TF_dict:
	#print "confirmed"
	#for x in confirmed_TF_dict[t]:
		#print x

#print "------------------------------"

#if t in tfdict:
	#print "results"
	#for item in tfdict[t]:
		#print item




	
			
"""
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/TF_GO_analysis.py')
"""

	
	
	
	


