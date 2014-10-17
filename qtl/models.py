from django.db import models

"""
Made by Jiao Long
"""

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
    strand = models.BooleanField(default=False,blank = True) # Django True: sense strand False: anti-sense strand or MySQL 1: sense strand 0: antisense strand
    
    #tair_accession = models.CharField(max_length=30) #Locus:2025361
    #gene__type = models.CharField(max_length=15) # protein_coding
    #associated_loci = models.TextField()
    #TEST FOR SYNCHRONIZATION
    
    def __unicode__(self):
        return self.locus_identifier
    
    
class Marker(models.Model):
    
    
    marker_name = models.CharField(max_length=15, primary_key=True) #PVV4.1
    marker_chromosome = models.IntegerField()#1
    marker_cm = models.DecimalField(max_digits = 3, decimal_places =1)#64.6
    marker_phys_pos = models.DecimalField(max_digits = 13, decimal_places =10)#0.008639
    #associated_locus = models.ForeignKey(Locus) #AT1G01480
    #marker_aliases = models.CharField(max_length=15) #PVV4
    #tair_accession = models.CharField(max_length=30) #GeneticMarker:1945638
    #marker_type = models.CharField(max_length=10) #CAPS
    #marker_length = models.CharField() #The data format of marker_length is like 1.000 (bp) 
    #is_PCR_marker = models.BooleanField() #true
    #special_condition = models.CharField()
    #chromosome = models.IntegerField() #1
    
    #map = models.ManyToManyField(Map) or add an instance in between to avoid of many-to-many relationship
        
    def __unicode__(self):
        return self.marker_name

class Experiment(models.Model):
    
    
    experiment_name=models.CharField(max_length=50,primary_key=True)
    
    def __unicode__(self):
        return self.experiment_name
    
class LOD(models.Model):
    
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gxe = models.BooleanField(default=False)
    locus_identifier = models.ForeignKey(Gene)
    marker_name = models.ForeignKey(Marker)
    
    def __unicode__(self):
        return self.LOD_score
    
class Parent(models.Model):
    
    parent_type = models.CharField(max_length=20)
    expression = models.DecimalField(max_digits = 25, decimal_places = 15)
    locus_identifier = models.ForeignKey(Gene)
    
    def __unicode__(self):
        return self.parent_type 

class RIL(models.Model):
    locus_identifier = models.ForeignKey(Gene)
    ril_name = models.CharField(max_length=20)
    ril_type = models.CharField(max_length=20)
    ril_exp = models.DecimalField(max_digits = 25, decimal_places = 15)
    def __unicode__(self):
        return self.ril_name 

class Metabolite(models.Model):
    metabolite_name = models.CharField(max_length=50,primary_key=True)
    def __unicode__(self):
        return self.metabolite_name 

class MParent(models.Model):
    parent_type = models.CharField(max_length=20)
    expression = models.DecimalField(max_digits = 25, decimal_places = 15)
    metabolite_name = models.ForeignKey(Metabolite)
    def __unicode__(self):
        return self.parent_type

class MRIL(models.Model):
    metabolite_name = models.ForeignKey(Metabolite)
    ril_name = models.CharField(max_length=20)
    ril_type = models.CharField(max_length=20)
    ril_exp = models.DecimalField(max_digits = 25, decimal_places = 15)
    def __unicode__(self):
        return self.ril_name

class MLOD(models.Model):
    experiment_name = models.ForeignKey(Experiment)
    LOD_score = models.DecimalField(max_digits = 12, decimal_places = 10)
    gxe = models.BooleanField(default=False)
    metabolite_name = models.ForeignKey(Metabolite)
    marker_name = models.ForeignKey(Marker)
    def __unicode__(self):
        return self.LOD_score


    



