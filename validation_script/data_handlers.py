import collections as col
import os

def read_data(fn):
	"""
	Eat and spit data from a file in a fast and memory-inefficient way
	"""
	with open(fn, 'r') as fo:
		data = fo.readlines()
	return data
	

def parse_AtReg_data(data):
	"""
	"""

	db = []
	for line in data:
		info = line.split('\t')
		
		nr = info[0].strip('"')
		gene_name = info[1].strip('"')
		TF = info[2].strip('"')
		family = info[3].strip('"')
		alt_name = info[4].strip('"')
		locus_id = info[5].strip('"')
		yesno = info[6].strip('"')
		direct = info[7].strip('"')
		confirmed = info[8].strip('"')
		biology = info[9].strip('"')
		activity = info[10].strip('"')
		article = info[11].strip('"')
		
		if len(TF.strip()) > 9:
			#print len(TF.strip())
			#print "unsplit: %s"%TF
			tf1 = TF[0:9]
			tf2 = TF[10:20]
			TF_split = [tf1, tf2]
			#print "split: %s"%TF_split
			for tf_n in TF_split:
				cell = [nr,gene_name,tf_n.upper(),family,alt_name,locus_id.upper(),yesno,direct,confirmed,biology,activity,article]
				db.append(cell)
		
		else:		
		
			cell = [nr,gene_name,TF.upper(),family,alt_name,locus_id.upper(),yesno,direct,confirmed,biology,activity,article]
			db.append(cell)
		
	
	return db
	

def parse_family_data(data):
	"""
	Retrieve a list with known Transciption Factors from a file
	The file was downloaded from AGRIS
	"""

	TF_list = []
	for line in data:
		info = line.split()
		TF_family = info[0]
		TF_locus = info[1].upper()
		TF_list.append([TF_family, TF_locus])
		
	
	return TF_list
	




def TFloc_pairs_AtRegNet(Dbase, select_fam):
	"""
	Retrieve a list of regulators which are linked to their regulated
	targets
	
	select_fam is a list of 1 or more TF family names
	"""
	
	all_families = [info[1] for info in Dbase]
	AtReg = []
	if len(select_fam) == 1:
		
		fam = select_fam[0]
		#Use all known regulator-target gene connections
		if fam == "all":
	
			for line in Dbase:

				TF = line[2].upper()
				loc = line[5].upper()
				
				if not [TF,loc] in AtReg:
					AtReg.append([TF,loc])
		
		elif fam not in all_families:
			print "The transcription family %s was not recognized"%fam
			
		else:
			#Use only the selected known regulator-target gene connection
			for line in Dbase:
				fam_name = line[1]
				
				if fam_name in select_fam and fam_name in all_families:
					
					TF = line[2].upper()
					loc = line[5].upper()
					
					if not [TF,loc] in AtReg:
						AtReg.append([TF,loc])			
			
		
	else:
		#Use all selected family names
		for line in Dbase:
			fam_name = line[1]
			
			if fam_name in select_fam and fam_name in all_families:
				
				TF = line[2].upper()
				loc = line[5].upper()
				
				if not [TF,loc] in AtReg:
					AtReg.append([TF,loc])
				
			elif fam_name in select_fam and fam_name not in all_families:
				print "The family name %s was not found in AtRegNet.txt"%fam
			
		
	return AtReg				


def make_AtReg_dict(data):
	"""
	Make a dict with Tf as key and target as value
	"""
	
	AtReg_dict = col.OrderedDict()
	for line in data:
		TF = line[2].upper()
		target = line[5].upper()
		
		if TF not in AtReg_dict:
			AtReg_dict[TF] = [target]
		else:
			temp_target = AtReg_dict[TF]
			temp_target.append(target)
			AtReg_dict[TF] = temp_target
		
	return AtReg_dict
		

	
	
def parse_verdict(verdict):
	"""
	Take the dictionary verdict and extract a list of regulators paired
	with their predicted targets
	
	verdict[(trait, go, region)] = [list of TFs]
	
	""" 
	
	regloc_list = []
	loc_list = []
	
	for key, value in verdict.iteritems():
		trait = key[0]
		#print "t",trait
		#print "k",key
		loc_list.append(trait)
		
		#populate the regloclist as such: [TF, locus]
		for TF in value:
			regloc_list.append([TF, trait])
			
	return regloc_list, set(loc_list)



def write_pairwise_validation(fn, data):
	"""
	Filename should contain information on dataset and cutoff
	"""
	
	with open(fn, 'w') as fo:
		for pair in data:
			fo.write("%s %s"%(pair[0], pair[1]))
			fo.write("\n")
	
	
def count_values(adict):
	cnts = 0
	for k, v in adict.iteritems():
		cnts += len(v)
	return cnts	
	

def create_dataset(TG_TF):
	"""
	"""
	tt_dict = col.OrderedDict()
	for tt in TG_TF:
		TG, TF = tt
		if TG not in tt_dict:
			tt_dict[TG] = [TF]
		else:
			tempTFlist = tt_dict[TG]
			tempTFlist.append(TF)
			tt_dict[TG] = tempTFlist
	
	return tt_dict
	
def make_folder(fn):
	"""
	"""
	if not os.path.exists(fn):
		os.mkdir(fn)
	
	
	
	
		

