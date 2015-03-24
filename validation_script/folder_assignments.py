#------------------------------------------------
#Set folder labels
#------------------------------------------------
#Main
main_folder = "/home/wouter/xenv_dj1.6/QTL"
script_folder = "validation_script"
result_folder = "validation"
mr_folder =  "%s/%s"%(main_folder, result_folder)

#Raw data
raw_folder = "raw_data"
filename_atreg = "%s/%s/AtRegNet.txt"%(mr_folder, raw_folder)
filename_fam = "%s/%s/families_data.tbl"%(mr_folder, raw_folder)

#Trait lists
trait_folder = "trait_lists"

#Enrichment results
enriched_folder = "enrichment_results"

#Validation results
validate_folder = "validation_results"

#TF analysis
TF_folder = "TF_analysis"
filename_TF = "%s/%s"%(mr_folder, TF_folder)#%(dataset)

#Scattered
scatplot_folder = "scatter_plots"
plots = "plots"

#Pairwise columns for TF and target genes -> validated results in talble form
pvfolder = "pairwise_validation"

gfolder = "gene_lists"

numfolder = "validation_numerics"

full_trait_info = "full_trait_info"



