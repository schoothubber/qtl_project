"""
Retrieve all traits with an eQTL for all datasets and various cutoffs
"""

import folder_assignments as fa
from data_handlers import read_data


def get_trait_with_genelist(dataset, cutoff, chromosome):
	"""
	Extract data from several files
	The files are arranged per dataset and cutoff
	
	The data consists of traits with a genelist from their eQTL regions
	"""
	
	subfolder = "genelist_%s"%dataset
	
	trait_list = []
			
	#Per trait extract the eQTL genelist
	for chrom in chromosome:

		filelocation = "%s/%s/%s/genelist_%s_co%s_chr%s.txt"%(
															fa.mr_folder, 
															fa.gfolder, 
															subfolder, dataset, 
															cutoff, chrom
															)


		#These files contain the genelists from the eQTLs of every trait
		with open(filelocation, 'r') as fo:
			for line in fo:
				
				if line.startswith('trait'):
					trait = line[7:16]
					trait_list.append(trait)
		
	return trait_list


def write_trait_to_file(fn, traits):
	"""
	"""
	print "writing %s to file"%fn
	with open(fn, 'w') as fo:
		for trait in traits:
			fo.write(trait)
			fo.write('\n')



def main():
	"""
		
	"""
		
	#Datasets
	exp_list = ['Ligterink_2014', 'Ligterink_2014_gxe', 'Snoek_2012','Keurentjes_2007']

	#Variables
	chromosome = [1,2,3,4,5]
	cutoff_list = [6.7, 4.3, 3]
	
	reduce_traits = False
	reduce_traits_even_more = True
	
	for dataset in exp_list:
		
		storage_folder = "%s/%s/genelist_%s"%(fa.mr_folder, fa.gfolder, dataset)

		for cutoff in cutoff_list:
			
			if reduce_traits:
				print "retrieving traits for %s %s"%(dataset, cutoff)
				traits = get_trait_with_genelist(dataset, cutoff, chromosome)
				
				fileloc = "%s/%s/reduced_traitlist_%s_co%s.txt"%(
																fa.mr_folder, 
																fa.trait_folder, 
																dataset, cutoff
																)
				
				print "writing file to %s"%fileloc												
				write_trait_to_file(fileloc, traits)
			
			if reduce_traits_even_more:
				#emr = even more reduced
				emr_traits = []
				fname = "%s/genelist_%s_co%s.txt"%(
									storage_folder, 
									dataset, cutoff
									)
				traitdata = read_data(fname)
				for line in traitdata:
					if line.startswith("trait:"):
						trait = line[7:16].strip()
						emr_traits.append(trait)
						
				#make a new emr_traitfilename
				#and store the traits
				
				emr_fileloc = "%s/%s/emr_traitlist_%s_co%s.txt"%(
																fa.mr_folder, 
																fa.trait_folder, 
																dataset, cutoff
																)
				print "writing file to %s"%emr_fileloc												
				write_trait_to_file(emr_fileloc, emr_traits)
	
	
	
	
main()
