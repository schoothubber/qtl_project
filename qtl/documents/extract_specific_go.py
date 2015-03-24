import csv
import string


def read_once_csv(filename):
	"""
	Read a csv file
	make a set of the first row
	"""
	
	alist = []
	
	with open(filename, 'rb') as csvobject:
		reader = csv.reader(csvobject)
		for row in reader:
			
			alist.append(row[0])
			
	aset = set(alist)
	
	return aset






def read_once_text(filename):
	
	fileobject = open(filename, 'r')
	data = fileobject.readlines()
	fileobject.close()
	
	
	return data




def read_go_data(data):
	"""

	"""
	
	go_list = []
	
	for gene_go in data:

		s = string.split(gene_go)

		if go_list.count(s[1]) == 0:
		
			go = s[1]
			go_list.append(go)
			
	go_set = set(go_list)

	return go_set


def main():
	
	AT_textfile = 'ArabiGOannotations.out'
	file_all = 'go_info.csv'
	
	data = read_once_text(AT_textfile)
	
	set_all = read_once_csv(file_all)
	set_reduced = read_go_data(data)
	
	newset = set_all & set_reduced
	
	print 'all = %d' %len(set_all)
	print 'reduced = %d' %len(set_reduced)
	print 'new = %d' %len(newset)
	
	



if __name__ == "__main__":
	
	main()
