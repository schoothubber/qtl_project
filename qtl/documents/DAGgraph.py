import pygraphviz as pgv


def draw_tree():
	G = pgv.Agraph(name="GO tree")
	edgeset = set()
	terms = []#list of go terms -> parents and child relations
	
	for term in terms:
		if draw_parents:
			edgeset.update(rec.get_all_parent_edges())
		if draw_children:
			edgeset.update(term.get_all_child_edges())

