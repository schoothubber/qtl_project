"""
Get eQTL size data
and plot it in a graph
"""

import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import folder_assignments as fa
from data_handlers import read_data



def process_data(data):
	"""
	"""
	sizes = []
	for line in data:
		if line.startswith("eQTL_size:"):
			size = int(line[11:].strip())
			sizes.append(size)
	return sizes
	

def data_plotter(fn, darray):
	"""
	"""
	fig = Figure(figsize = (10,10), dpi = 100)
	axis = fig.add_subplot(1,1,1)
	#axis.set_title(title)
	axis.set_xlabel("eQTL")
	axis.set_ylabel("size")
	axis.grid(False)
	axis.autoscale(enable = True)	
	
	axis.set_xlim(0, 4000)
	axis.set_ylim(0, 2000)
	#x = sorted(darray, reverse=True)
	x = darray
               
	N = len(x)

	ind = np.arange(N)  # the x locations for the groups
	width = 0.35       # the width of the bars
	
	axis.bar(ind*width, x, width, color = 'b')
	
	canvas = FigureCanvas(fig)
	print "Saving file..."
	canvas.print_figure(fn)	


def scatter_plotter(fn, data1, data2):
	"""

	"""
	fig = Figure(figsize = (10,10), dpi = 100)
	axis = fig.add_subplot(1,1,1)
	#axis.set_title(title)
	axis.set_xlabel("eQTL")
	axis.set_ylabel("size")
	axis.grid(False)
	axis.autoscale(enable = True)
	
	if len(data1) > len(data2):
		limit = len(data1)
	if len(data1) < len(data2):
		limit = len(data2)		

	axis.set_xlim(0, limit)
	axis.set_ylim(0, 2000)
	
	x = []
	y = []
	
	for nr, size in enumerate(data1):
		x.append(nr)
		y.append(size)
		
	axis.scatter(x,y, s=0.5, color='r', alpha=0.5)
	
	x = []
	y = []
	
	for nr, size in enumerate(data2):
		x.append(nr)
		y.append(size)
		
	axis.scatter(x,y, s=0.5, color='b', alpha=0.5)
	axis.legend( ("Ligterink_2014", "Keurentjes_2007"), loc = "upper right")
	
	canvas = FigureCanvas(fig)
	print "Saving file..."
	canvas.print_figure(fn)
	

def main():
	"""
	"""
	#options:
	bar_is_true = False
	scat_is_true = True

	plot_folder = "%s/%s/"%(fa.mr_folder, fa.plots)
	
	exp_list = ['Ligterink_2014', 'Keurentjes_2007']
	#exp_list = ['Keurentjes_2007']
	cutoff_list = [3]
	
	data_array = []
	
	for dataset in exp_list:
		
		for cutoff in cutoff_list:

			storage_folder = "%s/%s/eQTLsize_%s"%(fa.mr_folder, fa.gfolder, dataset)
			fname = "%s/eQTLsize_%s_co%s.txt"%(storage_folder, dataset, cutoff)

			raw_data = read_data(fname)
			data = process_data(raw_data)
			
			if bar_is_true:
				plot_fn = "%s/eQTL_plot_2%s_co%s.png"%(plot_folder, dataset, cutoff)
				data_plotter(plot_fn, data)
				
			if scat_is_true:
				data_array.append(data)


	data1 = data_array[0]
	data2 = data_array[1]
	
	scat_fn = "%s/eQTL_scatplot%s_co%s.png"%(plot_folder, dataset, cutoff)
	scatter_plotter(scat_fn, data1, data2)
	
	
main()








