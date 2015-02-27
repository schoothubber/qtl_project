import sys
import os
import networkx as nx
import matplotlib.pyplot as plt


from sortedcontainers import SortedDict

#The following lines enable this script to contact django
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')

def read_data(fn):
	"""
	Eat and spit data from a file in a fast and memory-inefficient way
	"""
	with open(fn, 'r') as fo:
		data = fo.readlines()
	return data
	

def parse_data(data):
	"""
	"""

	db = []
	for line in data:
		info = line.split('\t')
		
		nr = info[0].strip('"')
		binding = info[1].strip('"')
		TF = info[2].strip('"')
		family = info[3].strip('"')
		gene_name = info[4].strip('"')
		locus_id = info[5].strip('"')
		yesno = info[6].strip('"')
		direct = info[7].strip('"')
		confirmed = info[8].strip('"')
		biology = info[9].strip('"')
		activity = info[10].strip('"')
		article = info[11].strip('"')
		
		cell = [nr,binding,TF,family,gene_name,locus_id,yesno,direct,confirmed,biology,activity,article]
		db.append(cell)
		
		#articles.append(article)
		
		#if TF not in tf_dict:
			#tf_dict[TF] = [gene]
		#else:
			#tempgenes = tf_dict[TF]
			#tempgenes.append(gene)
			#tf_dict[TF] = tempgenes
	
	#return articles, tf_dict
	
	return db
	
def write_data(fn, dat):
	"""
	"""
	with open(fn, 'w') as fo:
		for info in dat:
			if info[1] == "AG":
			
				for item in info:
					fo.write("%s\t "%item)
				fo.write("\n\n")
				
				
def filter_data(dat, var):
	"""
	Select a range of locus ids based on their regulator
	the var variable will contain a string that indicates which
	regulator
	
	The nodes should be associated with an integer ranging from 0 to len(nodes)-1
	"""
	filtered = []
	for info in dat:
		if info[1] == var:
			filtered.append(info)
			
	edges = [(item[2], item[4]) for item in filtered]
	nodes = [item[4] for item in filtered]
	
	return edges, nodes

def parse_results(dat):
	"""
	"""
	results = []
	for line in dat:
		if line.startswith("trait:"):
			trait = line[7:16]
		if line.startswith("location:"):
			loc = line[10:13]
		if line.startswith("TF:"):
			TF = line[3:].split()
			
			for factor in TF:
				locus_factor_pair = [factor, trait]
				results.append(locus_factor_pair)
				
	return results


def hhh():
	"""
		Enter this command to open the django shell 
		python manage.py shell

		And enter this command to run this script	
		execfile('qtl/compare_with_literature.py')
	"""



	filein1 = 'raw_data/AtRegNet.txt'
	fileout = 'enrichment_results/AtRegNet_results.txt'

	AtReg_data = read_data(filein1)
	Dbase = parse_data(AtReg_data)
	#write_data(fileout, Dbase)


	filein2 = 'enrichment_results/TF_location_Ligterink_2014.txt'
	results_data = read_data(filein2)
	parse_res = parse_results(results_data)




	#----------------------------------------------------------------------
	#draw a network!
	#----------------------------------------------------------------------

	print "Initializing a graph, standby..."
	G = nx.Graph()


	print "Retrieving nodes, edges, labels and positions, standby..."

	edges, nodes = filter_data(Dbase, "AG")
	e_nodes = enumerate(nodes)

	#The labels 
	labels = {}
	for i in range(0, len(nodes)):
		labels[i] = nodes[i]
		
	#inverted labels (inverted dictionary)
	inv_labels = {v:k for k,v in labels.iteritems()}

	e_nodes = []
	for node in nodes:
		e_nodes.append(inv_labels[node])

	e_edges = []
	for edge in edges:
		temp_edge1 = inv_labels[edge[0]]
		temp_edge2 = inv_labels[edge[1]]
		print (temp_edge1, temp_edge2)
		e_edges.append((temp_edge1, temp_edge2))
		#G.add_edge(temp_edge1, temp_edge2)


	pos=nx.spectral_layout(G)
	print "----------------------------------------"
	print "The positions:"
	for k, v in pos.iteritems():
		print k, v
	print "----------------------------------------"
	print "The labels:"
	for k, v in labels.iteritems():
		print k, v
	print "----------------------------------------"


	# nodes
	nx.draw_networkx_nodes(G,pos, nodelist=e_nodes, node_color='b', node_size=100, alpha=0.8)

	# edges
	nx.draw_networkx_edges(G,pos, edgelist=e_edges, width=0.1, alpha=0.5, edge_color='r')
						   
	print "added %s edges to the graph"%len(edges)
	print "added %s nodes to the graph"%len(nodes)
	print "added %s labels to the graph"%len(labels)


	print "Commencing the buildup of the graph"


	nx.draw(G)
	nx.draw_networkx_labels(G, pos, labels,font_size=8)
	plt.show()

	print "Saving graph..."
	plt.savefig("enrichment_results/result_network.png")




#########################Scatter plots##################################
#-----------------------------------------------------------------------







