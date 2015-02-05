from django.db import models

class Gene(models.Model):
    
    locus_identifier = models.CharField(max_length=30,primary_key=True) #AT1G01480
    gene_model_name = models.CharField(max_length=30,blank = True)#AT1G01480.1
    gene_model_description = models.TextField(blank = True)
    gene_model_type = models.CharField(max_length = 40,blank = True)
    primary_gene_symbol = models.TextField(blank = True)
    all_gene_symbols = models.TextField(blank = True)
    chromosome = models.IntegerField(blank = True)
    start = models.IntegerField(blank = True)
    end = models.IntegerField(blank = True)
    orientation = models.BooleanField(blank = True) # Django True: sense strand False: anti-sense strand or MySQL 1: sense strand 0: antisense strand

    def __unicode__(self):
        return self.locus_identifier
    
class Marker(models.Model):
    
    marker_name = models.CharField(max_length=15, primary_key=True) #PVV4.1
    marker_chromosome = models.IntegerField()#1
    marker_cm = models.DecimalField(max_digits = 4, decimal_places =1)#101.6 
    marker_phys_pos = models.DecimalField(max_digits = 13, decimal_places =10)#0.008639

    def __unicode__(self):
        return self.marker_name

class Experiment(models.Model):

    experiment_name=models.CharField(max_length=50,primary_key=True)
    
    def __unicode__(self):
        return self.experiment_name
    
class ExperimentMarker(models.Model):
    experiment_name = models.ForeignKey(Experiment)
    marker_name = models.CharField(max_length=15)
    
class LOD(models.Model):
    
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gxe = models.BooleanField()
    locus_identifier = models.ForeignKey(Gene)
    marker_name = models.ForeignKey(Marker)
     
    def __unicode__(self):
        return self.LOD_score
    
class Parent(models.Model):
    
    parent_type = models.CharField(max_length=20)
    expression = models.DecimalField(max_digits = 25, decimal_places = 15)
    locus_identifier = models.ForeignKey(Gene)
    experiment_name = models.CharField(max_length=40,blank = True)
    
    def __unicode__(self):
        return self.parent_type 

class RIL(models.Model):
    locus_identifier = models.ForeignKey(Gene)
    ril_name = models.CharField(max_length=20)
    ril_type = models.CharField(max_length=20,blank = True)
    ril_exp = models.DecimalField(max_digits = 25, decimal_places = 15)
    experiment_name = models.CharField(max_length=40,blank = True)
    def __unicode__(self):
        return self.ril_name 


class Genotype(models.Model):
    ril_name = models.ForeignKey(Marker)
    ril_name = models.CharField(max_length=20)
    genotype = models.CharField(max_length=5,blank = True)# there might be some RIL populations missing genotype information. 
    experiment_name = models.CharField(max_length=40,blank = True)


class Metabolite(models.Model):
    metabolite_name = models.CharField(max_length=50,primary_key=True)
    def __unicode__(self):
        return self.metabolite_name 

class MParent(models.Model):
    parent_type = models.CharField(max_length=20)
    expression = models.DecimalField(max_digits = 25, decimal_places = 15)
    metabolite_name = models.ForeignKey(Metabolite)
    experiment_name = models.CharField(max_length=40,blank=True)
    def __unicode__(self):
        return self.parent_type

class MRIL(models.Model):
    metabolite_name = models.ForeignKey(Metabolite)
    ril_name = models.CharField(max_length=20)
    ril_type = models.CharField(max_length=20,blank = True)
    ril_exp = models.DecimalField(max_digits = 25, decimal_places = 15)
    experiment_name = models.CharField(max_length=40,blank = True)
    def __unicode__(self):
        return self.ril_name

class MLOD(models.Model):
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gxe = models.BooleanField()
    metabolite_name = models.ForeignKey(Metabolite)
    marker_name = models.ForeignKey(Marker)
    def __unicode__(self):
        return self.LOD_score
