import os
import sys

from django.test import TestCase
from django.utils import unittest

from qtl.models import Gene,Marker,LOD, Experiment

sys.path.append('home/env_dj1.6/QTL')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'QTL.settings')




class SimpleTest(TestCase):
	
    def test_basic_addition(self):
        """
        Tests that 1 + 1 always equals 2.
        """
        self.assertEqual(1 + 1, 2)
        self.assertEqual(45 / 9, 5)
        self.assertEqual(3 * 3, 9)

class MarkerLODscoreTests(TestCase):
	
	def test_individual_trait(self):
		"""
		Check if the url returns an HttpResponse
		And if the keys are passed along in the context
		
		Note that checking whether certain keys are present in the context
		of a response will always fail if there is no database connection
		"""
		Gene_1 = Gene.objects.create(
			locus_identifier = 'AT1G01010',
			gene_model_name = 'AT1G01010.1',
			gene_model_description = 'NAC domain containing protein 1 (NAC001); FUNCTIONS IN: sequence-specific DNA binding transcription factor activity; INVOLVED IN: multicellular organismal development, regulation of transcription; LOCATED IN: cellular_component unknown; EXPRESSED IN: 7 plant structures; EXPRESSED DURING: 4 anthesis, C globular stage, petal differentiation and expansion stage; CONTAINS InterPro DOMAIN/s: No apical meristem (NAM) protein (InterPro:IPR003441); BEST Arabidopsis thaliana protein match is: NAC domain containing protein 69 (TAIR:AT4G01550.1); Has 2503 Blast hits to 2496 proteins in 69 species: Archae - 0; Bacteria - 0; Metazoa - 0; Fungi - 0; Plants - 2502; Viruses - 0; Other Eukaryotes - 1 (source: NCBI BLink).',
			gene_model_type = 'protein_coding',
			primary_gene_symbol = 'NAC DOMAIN CONTAINING PROTEIN 1 (NAC001)',
			all_gene_symbols = 'ARABIDOPSIS NAC DOMAIN CONTAINING PROTEIN 1 (ANAC001);NAC DOMAIN CONTAINING PROTEIN 1 (NAC001)',
			chromosome = 1,
			start = 3631 ,
			end = 5899,
			strand = 1,
			)
			
		Marker_1 = Marker.objects.create(
			marker_name = 'ATHCHIB2',
			marker_chromosome = 3,
			marker_cm = 6.6,
			marker_phys_pos = 3.9633910000,
			)
			
		if Experiment.objects.filter(experiment_name = 'lod_test.txt').exists():
			print 'yes'
		else:
			print 'No'
			
		LOD_1 = LOD.objects.create(
			experiment_name = 'lod_test',
			LOD_score = -2.2960360000,
			locus_identifier = 'AT1G01010',
			marker_name = 'ATHCHIB2',
			gxe = 0,
			)
			
	
		
		response = self.client.get('/wouter/trait/')
		self.assertEqual(response.status_code, 200)
		self.assertTrue('link_for_trait' in response.context)
		self.assertTrue('Graph_url' in response.context)
		self.assertTrue('cutoff_name' in response.context)
		self.assertTrue('qtl_gogenes_dict' in response.context)
		self.assertTrue('go_genes_qtl_dict' in response.context)
		


































