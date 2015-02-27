import os
import collections as col




	
def write_verdict(verdict, dataset, cutoff, folder, subfolder):
	"""
	Write the results to a file
	"""	
	
	filename = "%s/%s/TF_results_%s_%s.txt"%(folder, subfolder, dataset, cutoff)
	
	if not os.path.exists(filename):

		with open(filename, 'w') as fo:
			for key, value in verdict.iteritems():
				#TGR stands for Trait, GO, Region
				fo.write("TGR: %s %s %s"%(key[0], key[1], key[2]))
				fo.write("\n")
				for tf in value:
					fo.write("%s "%tf)
				fo.write("\n\n")
	else:
		print "File %s already exists!"%filename
	



def verify(TF_set, filelocation):
	"""
	Enter this command to open the django shell 
	python manage.py shell

	And enter this command to run this script	
	execfile('qtl/verify_enrichment_results.py')
	"""
	#Get a list of known Transcription Factors
	#This list is compared to the enrichment results per trait
	#Any overlapping results will be recorded

	verify_dict = col.OrderedDict()
				
	with open(filelocation, 'r') as fo:
		for line in fo:
			
			if line.startswith('trait'):
				trait = line[7:16]
				
			if line.startswith('nr of eQTLS:'):
				eQTL = int(line[13:14])
				
			if line.startswith('GO:'):
				go = line.strip()
				
			if line[0].isdigit():
				region = line.strip()
				
			if line.startswith('AT'):
				enrichment = line.split()
				enrichment_set = set(enrichment)
				verification = enrichment_set & TF_set

				if verification:
					temp_tup = (trait, go, region)
					verify_dict[temp_tup] = list(verification)
	
						
	return verify_dict



					
#write_verdict(verify_dict, dataset, cutoff, fol, subfol)

					


"""
Enter this command to open the django shell 
python manage.py shell

And enter this command to run this script	
execfile('qtl/verify_enrichment_results.py')
"""
#folder = "enrichment_results"
#subfolder = "TF_analysis"
							
#verify(folder, subfolder)
				
	


				
