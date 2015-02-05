# qtl_project
Contains the necessary files for a Django app called 'qtl'

qtl/views.py is the controller for all app functions <br>
it imports functions from:<br>
genelist_from_eqtl.py (extracts a list of genes from the database based on eQTL cutoff values)<br>
go_enrichment.py (performs the functional annotation on the gene list)<br>
prepare_for_display.py (manipulates data for transfer to html sites)<br>
<br>
The templates are obviously stored in the folder 'templates'<br>
<br>
qtl/models.py is Django's template database design<br>

enrichment_results.py stores the said results based on cutoff and dataset.
