"""
Collect the Recalls, Precisions and F1 scores
for all datasets/cutoffs/eQTLthresholds
"""


import folder_assignments as fa
from data_handlers import read_data




def main():
	"""
	"""
	
	#-----------------------------------------------------------
	exp_list = ['Ligterink_2014','Ligterink_2014_gxe','Snoek_2012','Keurentjes_2007']
	cutoff_list = [3, 4.3, 6.7]
	eQTL_threshold = [0,1,2,3]
	chromo = [1,2,3,4,5]
	#-----------------------------------------------------------
	outputpath = "%s/%s/summary_validation_numerics.txt"%(
										 fa.mr_folder, fa.numfolder
											)	

	for eQTL in eQTL_threshold:
		
		for dataset in exp_list:
			
			substorage = "%s/%s/%s"%(fa.mr_folder, fa.numfolder, dataset)	
				
			for cutoff in cutoff_list:
			


				inputpath = "%s/permutate_eqtl_%s_%s_co%s.txt"%(
													 substorage, eQTL, 
													 dataset, cutoff
													)
				data = read_data(inputpath)
				
				for line in data:
					if line.startswith("lower_F1:"):
						lower_F1 = int(line[9:].strip())
					if line.startswith("higher_F1:"):
						higher_F1 = int(line[10:].strip())
						
					if line.startswith("lower_recall:"):
						lower_recall = int(line[13:].strip())
					if line.startswith("higher_recall:"):
						higher_recall = int(line[14:].strip())

					if line.startswith("lower_precision:"):
						lower_precision = int(line[16:].strip())
					if line.startswith("higher_precision:"):
						higher_precision = int(line[17:].strip())
							
				try:
					with open(outputpath, 'a') as fo:
						fo.write("dataset %s\n"%dataset)
						fo.write("cutoff %s\n"%cutoff)
						fo.write("eQTLs %s\n"%eQTL)
						fo.write("lower_F1 %s\n"%lower_F1)
						fo.write("higher_F1 %s\n"%higher_F1)
						fo.write("lower_recall %s\n"%lower_recall)
						fo.write("higher_recall %s\n"%higher_recall)
						fo.write("lower_precision %s\n"%lower_precision)
						fo.write("higher_precision %s\n"%higher_precision)
						fo.write("\n")
				except:
					pass


main()




