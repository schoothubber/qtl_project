


def parse_file(filename):
	"""
	Read an OBO file from the Gene Ontology Consortium.
	Downloaded from "http://geneontology.org/page/go-slim-and-subset-guide"
	"""
	
	
	with open(filename, 'r') as fileobject:
		
		data = fileobject.readlines()
		
	return data


def get_go_definition(data):
	"""
	Extract the Go descriptions from the file
	"""
	
	go_def_dict = {}
	i = 0
	for line in data:
		
		
		if line.startswith('id: GO'):
			goterm = line[4:].strip()
			i += 1
			print goterm
			print i
			
		#elif line.startswith('def:'):
			#definition = line[5:].strip()
			
		#go_def_dict[goterm] = definition
		
	#return go_def_dict





if __name__ == "__main__":
	
	filename = "goslim_plant.obo"
	filename2 = "goslim_generic.obo"
	
	
	d = parse_file(filename2)
	
	get_go_definition(d)
	#adict = get_go_definition(d)
	
	#for k, v in adict.iteritems():
		#print "For go term %s I found this definition %s" %(k, v)
	
	
	
	
	
	
	
