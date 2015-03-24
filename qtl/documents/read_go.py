import csv

def read_go_file(filename):
	
	go_name_tup_list = []
	
	go = None
	name = None
	namespace = None
	definition = None
	
	with open(filename, 'r') as fo:
		
		for line in fo:
		
			if line.startswith('id:'):
				
				go = line[3:].strip()
		
			elif line.startswith('name:'):
				
				name = line[5:].strip()
				
			elif line.startswith('namespace:'):
				
				namespace =  line[10:].strip()
				
			elif line.startswith('def:'):
				
				definition = line[4:].strip()
				
			elif line.startswith('[Term]'):
					
				go_name_tup = (go, name, namespace, definition)
				go_name_tup_list.append(go_name_tup)
		
	return go_name_tup_list
	
	
	
def split_on_namespace(tup_list):
	
	#define namespaces
	tup_list_bio = []
	tup_list_mol = []
	tup_list_cel = []
	
	
	for tup in tup_list:
		
		if tup[2] == 'biological_process':
			tup_list_bio.append(tup)
		
		if tup[2] == 'molecular_function':
			tup_list_mol.append(tup)
		
		if tup[2] == 'cellular_component':
			tup_list_cel.append(tup)
			
	return tup_list_bio, tup_list_mol, tup_list_cel
		


def write_list_to_csv(array, filename):
	
	with open(filename, 'wb') as csvfile:
		
		for alist in array:
		
			writer = csv.writer(csvfile)
			writer.writerow(alist)



def main():
	"""
	Read all go terms from the parsed obo file
	
	Extract and write desired data from text to csv file
	
	data(goterm, name, namespace, definition)
	"""
	fn_bio = 'go_bio_process.csv'
	fn_mol = 'go_mol_function.csv'
	fn_cel = 'go_cel_component.csv'
	fn_go = 'go_names.txt'
	fn_out_csv = 'go_info.csv'
	

	#make data
	tup_list = read_go_file(fn_go)
	
	#separate data
	bio, mol, cel = split_on_namespace(tup_list)
	
	#write data
	write_list_to_csv(bio, fn_bio)
	write_list_to_csv(mol, fn_mol)
	write_list_to_csv(cel, fn_cel)
	
	
	#write_list_to_csv(gn_list, fn_out_csv)




if __name__ == "__main__":
	
	main()
	

	
	
	
	
	
