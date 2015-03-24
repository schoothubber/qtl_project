import sys
import os

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from qtl.models import Metabolite, Marker,MLOD, Experiment,ExperimentMarker


##Set the path so that the functions in this file can be imported by views.py
sys.path.append('/home/wouter/xenv_dj1.6/QTL')
#sys.path.append('/mnt/geninf15/prog/www/django/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')


#Main
main_folder = "/home/wouter/xenv_dj1.6/QTL"
script_folder = "validation_script"
result_folder = "validation"
meta_plot_folder = "metabolite_plots"
mr_folder =  "%s/%s"%(main_folder, result_folder)


#class Marker(models.Model):
    
    #marker_name = models.CharField(max_length=15, primary_key=True) #PVV4.1
    #marker_chromosome = models.IntegerField()#1
    #marker_cm = models.DecimalField(max_digits = 4, decimal_places =1)#101.6 
    #marker_phys_pos = models.DecimalField(max_digits = 13, decimal_places =10)#0.008639

    #def __unicode__(self):
        #return self.marker_name

#class Experiment(models.Model):

    #experiment_name=models.CharField(max_length=50,primary_key=True)
    
    #def __unicode__(self):
        #return self.experiment_name
    
#class ExperimentMarker(models.Model):
    #experiment_name = models.ForeignKey(Experiment)
    #marker_name = models.CharField(max_length=15)
    
#class Metabolite(models.Model):
    #metabolite_name = models.CharField(max_length=50,primary_key=True)
    #def __unicode__(self):
        #return self.metabolite_name 

#class MParent(models.Model):
    #parent_type = models.CharField(max_length=20)
    #expression = models.DecimalField(max_digits = 25, decimal_places = 15)
    #metabolite_name = models.ForeignKey(Metabolite)
    #experiment_name = models.CharField(max_length=40,blank=True)
    #def __unicode__(self):
        #return self.parent_type

#class MRIL(models.Model):
    #metabolite_name = models.ForeignKey(Metabolite)
    #ril_name = models.CharField(max_length=20)
    #ril_type = models.CharField(max_length=20,blank = True)
    #ril_exp = models.DecimalField(max_digits = 25, decimal_places = 15)
    #experiment_name = models.CharField(max_length=40,blank = True)
    #def __unicode__(self):
        #return self.ril_name

#class MLOD(models.Model):
    #experiment_name = models.ForeignKey(Experiment)
    #LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    #gxe = models.BooleanField()
    #metabolite_name = models.ForeignKey(Metabolite)
    #marker_name = models.ForeignKey(Marker)
    #def __unicode__(self):
        #return self.LOD_score




"""
python manage.py shell
execfile('qtl/metabolite_mQTL.py')
"""

def get_metabolites():
	"""
	"""
	met_list = Metabolite.objects.values_list('metabolite_name', flat = True)
	
	return met_list
	

def get_markers(metabo, experiment, gxe_boolean):
	"""
	"""
	
	markers_full = ExperimentMarker.objects.filter(experiment_name = experiment)
	markers = markers_full.values_list('marker_name', flat = True)
	
	#for m in markers:
		#print m
		
	marker_logp_tuple = ()
	tuple_list = []
	 
	for i in range(0, len(markers)):
		
		try:
			logp1 = MLOD.objects.get(
								metabolite_name = metabo, marker_name = markers[i], 
								gxe = gxe_boolean, experiment_name = experiment
								)
								
			logp2 = logp1.LOD_score
			
			logp3 = float(logp2)
			
			marker_logp_tuple = (i, markers[i].encode("utf-8"), logp3)
			tuple_list.append(marker_logp_tuple)
	        
		except MLOD.DoesNotExist:
			pass
	
	return tuple_list





def draw_plot(meta, data1, data2, fn):
	"""
	
	"""
	

	fig = Figure(figsize = (10,10))
	axis = fig.add_subplot(1,1,1)
	axis.set_title(meta)
	axis.set_xlabel("MLOD")
	axis.set_ylabel("Markers")
	axis.grid(False)
	axis.autoscale(enable = True)	
	axis.axhline(0, color = 'k')	
	
	x = [info[0] for info in data1]
	y = [info[2] for info in data1]	
	axis.plot(x, y, color='r', label = 'sse')
	
	x = [info[0] for info in data2]
	y = [info[2] for info in data2]	
	axis.plot(x, y, color='b', label = 'gxe')
	
	box = axis.get_position()
	axis.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	
	axis.legend(loc = 'center left', bbox_to_anchor=(1.0,0.5))
	
	canvas = FigureCanvas(fig)
	print "Saving file..."
	canvas.print_figure(fn)


########################################################################
########################################################################
########################################################################
"""
python manage.py shell
execfile('qtl/metabolite_mQTL.py')
"""

cutoff_list = [3, 4.3, 6.7]
exp_list = ['Ligterink_2014','Ligterink_2014_gxe','Keurentjes_2007','Snoek_2012']


metabs = get_metabolites()


#draw a plot for every metabolite in the Database
#each plot will contain the steady state and the environmental variation
for umet in metabs:
	metabo = umet.encode("utf-8")
	filename = "Graph_%s"%metabo
	
	experiment = exp_list[0]
	
	gxe_boolean = False
	plotfilelocation = "%s/%s/%s.png"%(mr_folder, meta_plot_folder, filename)

	tup1 = get_markers(metabo, experiment, gxe_boolean)

	gxe_boolean = True
	tup2 = get_markers(metabo, experiment, gxe_boolean)

	draw_plot(metabo, tup1, tup2, plotfilelocation)








