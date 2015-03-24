import collections as col

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import numpy as np



def get_all_loci(main, sub, chromosome):
	"""
	Retrieve full list of locus ids
	"""
	complete_locus_list = []
	chrom_length_list = []
	chrom_start_pos = 0
	
	for chrom in range(1,chromosome+1):
		traitfile = '%s/%s/traitlist_chr%s.txt'%(main, sub, chrom)
		with open(traitfile, 'r') as fo:
			data = fo.readlines()
		traitlist = [trait.strip() for trait in data]
		complete_locus_list+=traitlist
		
		#calculate the starting points of the chromosomes
		#On these position there will be a yellow vertical line
		chrom_length = len(traitlist) + chrom_start_pos
		chrom_start_pos = chrom_length
		chrom_length_list.append(chrom_length)
	
	locus_dict = col.OrderedDict()
	for index, locus in enumerate(sorted(complete_locus_list)):
		locus_dict[locus] = index
		
	return locus_dict, chrom_length_list

	
def get_loci_from_dataset(locus_list):
	"""
	Retrieve a list of loci that were mapped
	i.e.: with 1, 2, 3 or more eQTLs
	These will indicate values on the x-axis of a scatterplot
	"""
	
	#print "len loc %s"%len(locus_list)
	
	
	#create the x-axis values			
	locus_dict = col.OrderedDict()
	for index, locus in enumerate(sorted(locus_list)):
		locus_dict[locus] = index
		
	#Calculate chromosomal regions of the TF list
	#these will be incicated by yellow vertical lines
	chromlist = [gene[2] for gene in locus_list]
	count_chrom = col.Counter(chromlist)

	ctf_list = []
	chrom_length_list = []
	chrom_start_pos = 0
	for k, v in count_chrom.iteritems():
		ctf_list.append([k, v])
		
	sorted_ctf_list = sorted(ctf_list)
	#Calculate starting posistions
	for chrom, tally in sorted_ctf_list:
		chrom_start_pos += tally
		chrom_length_list.append(chrom_start_pos)
		
	return locus_dict, chrom_length_list


def create_reg_dict(reg_list):
	"""
	Get a list of lists containing all known regulators
	These will come on the Y-axis
	The data is taken from families_data.tbl
	Spit out a dictionary of index zero with an enumeration
	"""
	chrom_length_list = []
	chrom_start_pos = 0
	
	reg_list = sorted(reg_list)
	#Calculate chromosomal regions of the TF list
	#these will be incicated by yellow horizontal lines
	chromlist = [gene[2] for gene in reg_list]
	count_chrom = col.Counter(chromlist)
	#place these  counted values in a sorted list
	ctf_list = []
	for k, v in count_chrom.iteritems():
		ctf_list.append([k, v])
	sorted_ctf_list = sorted(ctf_list)
	#Calculate chromosomal positions
	for chrom, tally in sorted_ctf_list:
		chrom_start_pos += tally
		chrom_length_list.append(chrom_start_pos)
		
	#create a list of enumerated TFs
	complete_reg_set = set(reg_list)
	reg_list = list(complete_reg_set)
	#standardize the list
	#stan_reglist = [reg.upper() for reg in reg_list]
	#sort the list
	sort_reg_list = sorted(reg_list)
	
	
	reg_dict = col.OrderedDict()
	for index, reg in enumerate(sort_reg_list):
		#print "%s - %s"%(index,reg)
		reg_dict[reg] = index
		
	return reg_dict, chrom_length_list
	
	

def link_data_in_array(data, locus_dict, reg_dict, AtReg_dict):
	"""
	data contains the pairwise relations between TFs and Loci from
	the model.
	The array will have size len(locus_dict) x len(reg_dict)
	'0' indicates no connection
	'1' indicates a connection
	"""
	
	xaxis = len(reg_dict)
	yaxis = len(locus_dict)
	size = (yaxis, xaxis)
	
	#Create a matrix of zeros
	data_array = np.zeros(size)


	
	print "matrix size is",size	
	
	labels = []
	#Populate array with meaningfull data:
	for i in range(0, len(data)):
		loc = data[i][0]
		reg = data[i][1]
		
		#if there is an x and y value in the matrix for this paired
		#TF-gene relationship:
		if reg in reg_dict and loc in locus_dict:
			#print loc
			
			#AtReg_dict can be empty of filled. If its filled than the 
			#connections in the pairwise columns from the model will be
			#compared to the connections from AtRegNet
			if AtReg_dict:
		
				if reg in AtReg_dict and loc in AtReg_dict[reg]:
					
					index_reg = reg_dict[reg]
					index_loc = locus_dict[loc]
					data_array[index_reg, index_loc] = 1
					
				#Creation of Labels for the scatterplot
					loc_label_data = (index_loc, index_reg, reg, loc)
					if loc_label_data not in labels:
						labels.append(loc_label_data)
			
			#If AtRegNet is empty than the validation is less strict:				
			else:
				
				index_reg = reg_dict[reg]
				index_loc = locus_dict[loc]
				data_array[index_loc, index_reg] = 1				
			
			#Creation of Labels for the scatterplot
				loc_label_data = (index_loc, index_reg, reg, loc)
				if loc_label_data not in labels:
					labels.append(loc_label_data)
		
		
	print "size of array is %s by %s"%(len(data_array), len(data_array[0]))
	
	
	return data_array, labels
	

def draw_plot(darray, fn, title, chr_len_list, TF_chr_len_list, labels, lab_lev, col):
	"""
	
	"""
	
	fig = Figure(figsize = (10,10), dpi = 100)
	axis = fig.add_subplot(1,1,1)
	axis.set_title(title)
	axis.set_xlabel("Transcription Factor")
	axis.set_ylabel("Target Genes")
	axis.grid(False)
	axis.autoscale(enable = True)	
	
	axis.set_xlim(0, len(darray[0]))
	axis.set_ylim(0, len(darray))
		
	for chrom_len in chr_len_list:
		axis.axhline(chrom_len, color = 'y', alpha = 0.5)
		
	for TF_chrom_len in TF_chr_len_list:
		axis.axvline(TF_chrom_len, color = 'y', alpha = 0.5)
			
	x, y = populate_the_axis(darray)
	axis.scatter(x,y, s=0.5, color=col, alpha=0.5)
	
	if labels and lab_lev:
	#add labels (x, y, TFlabel, LOClabel)
		for x, y, TF, loc in labels:
			#print "x is %s and y is %s and TF is %s and loc is %s"%(x, y, TF, loc)
			
			axis.annotate(loc, xy=(x,y), xycoords='data', xytext=(0,45), alpha=0.8, textcoords='offset points', 
					size=6, va="center", rotation=45,
					#bbox=dict(boxstyle="round4", fc="w"),
					  arrowprops=dict(arrowstyle="-|>",
									  connectionstyle="arc3,rad=0.2",
									  relpos=(0., 0.),
									  fc="w")
									)
									
			axis.annotate(TF, xy=(x,y), xycoords='data', xytext=(0,-45), alpha=0.8, textcoords='offset points', 
			size=6, va="center", rotation=45,
			#bbox=dict(boxstyle="round4", fc="w"),
			  arrowprops=dict(arrowstyle="-|>",
							  connectionstyle="arc3,rad=0.2",
							  relpos=(0., 0.),
							  fc="w")
						)                        

	canvas = FigureCanvas(fig)
	print "Saving file..."
	canvas.print_figure(fn)


def draw_combi_plot(darray, rarray, fn, title, chr_len_list, TF_chr_len_list, labels, lab_lever):
	"""
	darray holds data for the model
	rarray holds data for the reference (AtRegNet.txt
	
	"""
	

	fig = Figure(figsize = (10,10), dpi = 100)
	axis = fig.add_subplot(1,1,1)
	axis.set_title(title)
	axis.set_xlabel("Target Genes")
	axis.set_ylabel("Transcription Factor")
	axis.grid(False)
	axis.autoscale(enable = True)
	
	axis.set_xlim(0, len(darray[0]))
	axis.set_ylim(0, len(darray))
	
	if chr_len_list:
		for chrom_len in chr_len_list:
			axis.axvline(chrom_len, color = 'y', alpha = 0.5)
			
	if TF_chr_len_list:
		for TF_chrom_len in TF_chr_len_list:
			axis.axhline(TF_chrom_len, color = 'y', alpha = 0.5)
			
	if len(darray) != 0:
		x, y = populate_the_axis(rarray)
		axis.scatter(x,y, s=0.5, color='r', alpha=0.8)

	if len(rarray) != 0:
		x, y = populate_the_axis(darray)
		axis.scatter(x,y, s=1, color='b', alpha = 1)
	
	if labels and lab_lever:
	#add labels (x, y, TFlabel, LOClabel)
		for x, y, TF, loc in labels:
			#print "x is %s and y is %s and TF is %s and loc is %s"%(x, y, TF, loc)
			
			axis.annotate(loc, xy=(x,y), xycoords='data', xytext=(0,45), alpha=0.8, textcoords='offset points', 
					size=6, va="center", rotation=45,
					#bbox=dict(boxstyle="round4", fc="w"),
					  arrowprops=dict(arrowstyle="-|>",
									  connectionstyle="arc3,rad=0.2",
									  relpos=(0., 0.),
									  fc="w")
									)
	
	
	canvas = FigureCanvas(fig)
	
	print "Saving file..."
	canvas.print_figure(fn)


def populate_the_axis(data):
	"""
	Take data and give back x and y coordinates
	"""	
	x = []
	y = []
	
	for ind_1, sublist in enumerate(data):
		#print "data: ",data[:10]
		for ind_2, ele in enumerate(sublist):
			#print "sub: ",sublist[:10]
			if ele == 1:
				x.append(ind_2)
				y.append(ind_1)
				
	return x, y



#def activate():
	#"""
	#Enter this command to open the django shell 
	#python manage.py shell

	#And enter this command to run this script	
	#execfile('qtl/vis_regulator.py')
	#"""

	#folder = "enrichment_results"
	#subfolder = "TF_analysis"

	#scatfn1 = "scatter_atreg2.png"
	#scatfn2 = "scatter_model2.png"


	##Data from all traits
	#loc_dict = get_all_loci(5)

	##And from all regulators
	#reg_dict = create_reg_dict(get_family_list())
	#print "REGGG: %s"%len(reg_dict)

	##Data from AtReg
	#AtReg = get_all_regulators()
	#data_array1 = link_data_in_array(AtReg, loc_dict, reg_dict)

	##Data from model:
	#dresults = get_model_results()
	#data_array2 = link_data_in_array(dresults, loc_dict, reg_dict)

	##make plot:
	#print "Saving scatter plots..."
	#draw_plot(data_array1, scatfn1)

	#draw_plot(data_array2, scatfn2)





















	

