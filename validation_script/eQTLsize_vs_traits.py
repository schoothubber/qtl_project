"""
plot the number of gene in the eQTL against the traits
"""

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import folder_assignments as fa
from data_handlers import get_genelist, read_data, get_trait_samplesize_data, get_info


			
			
def draw_plot(x, y, fn, title):
	"""
	
	"""
	
	fig = Figure(figsize = (10,10), dpi = 100)
	axis = fig.add_subplot(1,1,1)
	axis.set_title(title)
	axis.set_xlabel("Traits")
	axis.set_ylabel("eQTL size")
	axis.grid(False)
	axis.autoscale(enable = True)	
	
	#axis.set_xlim(0, len(x))
	#axis.set_ylim(0, 5000)

	#for TF_chrom_len in TF_chr_len_list:
		#axis.axvline(TF_chrom_len, color = 'y', alpha = 0.5)
			
	axis.bar(x,y,)
	#axis.set_color('b')
	                  

	canvas = FigureCanvas(fig)
	print "Saving file..."
	canvas.print_figure(fn)			
			

def draw_3d_plot(x, y, z, fn, title):
	"""
	"""
	fig = Figure(figsize = (10,10), dpi = 100)
	axis = fig.add_subplot(111, projection='3d')
	axis.set_title(title)
	axis.set_xlabel("Traits")
	axis.set_ylabel("nr of eQTLs")
	axis.set_zlabel("genes in eQTLs")
	axis.grid(False)
	axis.autoscale(enable = True)	
	
	axis.set_xlim(0, len(x))
	#axis.set_ylim(0, 5000)

	#for TF_chrom_len in TF_chr_len_list:
		#axis.axvline(TF_chrom_len, color = 'y', alpha = 0.5)
			
	axis.bar(x,y,z)
	#axis.set_color('b')
	                  

	canvas = FigureCanvas(fig)
	print "Saving file..."
	canvas.print_figure(fn)		


	
def write_plot(tlist, fn, title):
	"""
	"""
	
	with open(fn, 'w') as fo:
		fo.write("%s"%title)
		fo.write("\n")		
		fo.write("-----------------------")
		fo.write("\n")
		for g, t in tlist:
			fo.write("%s\t%s"%(t, g))
			fo.write("\n")
			



		



def main():
	"""
	"""

	#exp_list = ['Ligterink_2014', 'Ligterink_2014_gxe', 'Keurentjes_2007', 'Snoek_2012']
	#cutoff_list = [3, 4.3, 6.7]
	#-----------------------------------------------------------
	
	#exp_list = ['Ligterink_2014']
	#exp_list = ['Ligterink_2014_gxe']
	#exp_list = ['Keurentjes_2007']
	exp_list = ['Snoek_2012']
	
	cutoff_list = [3]
	#cutoff_list = [4.3]
	#cutoff_list = [6.7]
	
	chromo = [1,2,3,4,5]
	#-----------------------------------------------------------
	
	fileloc = "%s/%s/tt_te_combi.txt"%(fa.mr_folder, fa.numfolder)
	#print "Retrieving random sample sizes from %s"%fileloc
	szdata = read_data(fileloc)
	sample_size_dict = get_trait_samplesize_data(szdata)
	
	draw_trait_vs_eqtlsize = False
	draw_trait_vs_nreqtls = True
	
	for dataset in exp_list:
		
		for cutoff in cutoff_list:
			
			key = (dataset, cutoff)	
			if key in sample_size_dict:
				#sample_size_list = [trait, sample_size]
				sample_size_list = sample_size_dict[key]
				tt_genes = [item[0] for item in sample_size_list]
			
########################################################################			
			x = []
			y = []	
########################################################################
			

			
			#Print some specific traits based on size of genelist
			#for t, gl in trait_genelist_list:
				#if len(gl) < 124:
					#print t
########################################################################
			
			if draw_trait_vs_eqtlsize:
				
				#trait_genelist_list = [trait, genelist]
				trait_genelist_list = get_genelist(dataset, cutoff, chromo)
				
				tt_trait_genelist_list = [[t[0], t[1]] for t in trait_genelist_list if t[0] in tt_genes]
				t_gls_list = [[len(info[1]), info[0]] for info in tt_trait_genelist_list]
				sort_t_gls_list = sorted(t_gls_list, reverse=True)
				t = [info[1] for info in sort_t_gls_list]
				gls = [info[0] for info in sort_t_gls_list]
				
				for ind_t, trait in enumerate(t):
					x.append(ind_t)
				
				for eQTLsize in gls:
					y.append(eQTLsize)
			
				filename_plot = "%s/plots/tt_eQTLvsTrait_%s_co_%s.png"%(fa.mr_folder, dataset, cutoff)
				filename_text = "%s/plots/tt_eQTLvsTrait_%s_co_%s.txt"%(fa.mr_folder, dataset, cutoff)
				title = "eQTL size vs Traits for %s with cutoff %s"%(dataset, cutoff)
			
			
				write_plot(sort_t_gls_list, filename_text, title)
				draw_plot(x, y, filename_plot, title)
				
				tot = float(sum(gls))
				if len(gls) != 0:
					avg = tot/float(len(gls))
				
					print key
					print "Average eQTL size expressed in nr of genes: %s"%avg
				
				
			if draw_trait_vs_nreqtls:		
				tr_eqtl_list = get_info(dataset, cutoff, chromo, tt_genes)

				tt_trait_eqtls_list = [[t[0], t[1]] for t in tr_eqtl_list if t[0] in tt_genes and t[1]>2]
				t_eqtls_list = [[info[1], info[0]] for info in tt_trait_eqtls_list]
				sort_t_eqtls_list = sorted(t_eqtls_list, reverse=True)
				t = [info[1] for info in sort_t_eqtls_list ]
				eqtls = [info[0] for info in sort_t_eqtls_list ]

				for ind_t, trait in enumerate(t):
					x.append(ind_t)
				
				for eQTLsize in eqtls:
					y.append(eQTLsize)
			
				filename_plot = "%s/plots/tt_nrofeQTLsvsTrait_%s_co_%s.png"%(fa.mr_folder, dataset, cutoff)
				filename_text = "%s/plots/tt_nrofeQTLsvsTrait_%s_co_%s.txt"%(fa.mr_folder, dataset, cutoff)
				title = "nr of eQTLs vs Traits for %s with cutoff %s"%(dataset, cutoff)
			
			
				#write_plot(sort_t_gls_list, filename_text, title)
				draw_plot(x, y, filename_plot, title)				
			
			

			
			
			
			

main()
