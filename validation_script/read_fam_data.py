
from collections import Counter, OrderedDict

from data_handlers import read_data, parse_AtReg_data, parse_family_data, TFloc_pairs_AtRegNet, parse_verdict



#------------------------------------------------
#Set folder labels
#------------------------------------------------
#Main
main_folder = "/home/wouter/xenv_dj1.6/QTL"
script_folder = "validation_script"
result_folder = "validation"
mr_folder =  "%s/%s"%(main_folder, result_folder)

#Raw data
raw_folder = "raw_data"
filename_atreg = "%s/%s/AtRegNet.txt"%(mr_folder, raw_folder)
filename_fam = "%s/%s/families_data.tbl"%(mr_folder, raw_folder)

filename_counted_fams = "%s/counted_fams.txt"%mr_folder
filename_counted_regs = "%s/counted_regs.txt"%mr_folder

#----------------------------------------------------------------------
#Set all data files and lists
#----------------------------------------------------------------------
#read data from files
TF_family_data = read_data(filename_fam)
AtReg_data = read_data(filename_atreg)
#parse data from file objects
TF_fam_data = parse_family_data(TF_family_data)
AtRegNet_parse = parse_AtReg_data(AtReg_data)
#make single list of TFs

TF_fam_list = [info[1] for info in TF_fam_data]
#turn these lists into into sets
#AtRegNet_set = set(AtRegNet_list)
#TF_fam_set = set(TF_fam_list)


AtRegNet_list = [info[2] for info in AtRegNet_parse]

for TF in AtRegNet_list:
	
	if len(TF.strip()) > 9:
		print TF



def count_regs(fn, AtRegNet_parse):
	"""
	"""
	AtRegNet_list = [info[2] for info in AtRegNet_parse]

	
	AtReg_dict = OrderedDict()
	for line in AtRegNet_parse:
		nam = line[1]
		loc = line[2]
		if loc not in AtReg_dict:
			AtReg_dict[loc] = nam

			
	
	cregs = Counter(AtRegNet_list)
	cregs_list = []
	for k, v in cregs.iteritems():
		cregs_list.append([v, k])
			
	sorted_cregs_list = sorted(cregs_list, reverse = True)
	
	with open(fn, 'w') as fo:
		for cnt, loc in sorted_cregs_list:
			fo.write("%s\t\t"%loc)
			if loc in AtReg_dict:
				fo.write("%s\t\t"%AtReg_dict[loc])
			fo.write("%s"%cnt)
			fo.write("\n")




def famcounter(fn, TF_fam_data):
	"""
	"""
	TF_fam_list = [info[0] for info in TF_fam_data]
	
	cfam = Counter(TF_fam_list)
	cfam_list = []
	for k, v in cfam.iteritems():
		cfam_list.append([v, k])
		
	sorted_cfam_list = sorted(cfam_list, reverse = True)

	with open(fn, 'w') as fo:
		for x in sorted_cfam_list:

			fo.write("%s - %s "%(x[1], x[0]))
			fo.write("\n")


def parse_more_family_data(data):
	"""
	Retrieve a list with known Transciption Factors from a file
	The file was downloaded from AGRIS
	"""

	TF_list = []
	for line in data:
		info = line.split()
		
		TF_family = info[0]
		TF_locus = info[1].upper()
		TF_name = info[2]
		
		TF_list.append([TF_locus, TF_family, TF_name])
		
	return TF_list



def famlist_per_chromosome(fn, TF_fam_list, TF_info_list):
	"""
	"""
	TF_list = sorted(TF_fam_list)
	
	chromlist = [gene[2] for gene in TF_list]
	count_chrom = Counter(chromlist)
	ctf_list = []
	for k, v in count_chrom.iteritems():
		ctf_list.append([k, v])	
		
		
	sTF_info_list = sorted(TF_info_list)
		
	
	with open(fn, 'w') as fo:
		
		fo.write("%s"%len(TF_list))
		fo.write("\n")
		for chrom, nr in ctf_list:
			fo.write("chromosome %s has %s TFs"%(chrom, nr))
			fo.write("\n")
		fo.write("\n\n")
		fo.write("-------------------------------------\n")		
		for loc, fam, des in sTF_info_list:

			fo.write("%s\t%s\t\t\t%s"%(loc, fam, des))
			fo.write("\n")





#TFinfo_list = parse_more_family_data(TF_family_data)
#filename_fpc = "%s/TF_per_chromosome.txt"%mr_folder
#famlist_per_chromosome(filename_fpc, TF_fam_list, TFinfo_list)


#count_regs(filename_counted_regs, AtRegNet_parse)


