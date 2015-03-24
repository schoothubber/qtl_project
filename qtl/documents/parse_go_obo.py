
#parsing http://geneontology.org/ontology/go.obo


import urllib2


def read_page(url_adress):
	response = urllib2.urlopen(url_adress)
	data = response.read()
	return data
	
def save_data(data, filename):
	with open(filename, 'w') as fileobject:
		
		for line in data:
			fileobject.write(line)





if __name__ == "__main__":
	
	url = 'http://geneontology.org/ontology/go.obo'
	
	fn = 'go_names.txt'
	
	
	GO_data = read_page(url)
	
	save_data(GO_data, fn)
