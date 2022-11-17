# ht_prs_preg
R code for the publication: Associations of polygenic risk scores for preeclampsia and blood pressure
with hypertensive disorders of pregnancy
https://doi.org/10.1097/HJH.0000000000003336

PRS values were calculuted using PRS-CS pipeline with default settings: https://github.com/getian107/PRScs. Calculations were performed at the FinnGen PRS pipeline: https://github.com/FINNGEN/CS-PRS-pipeline

```
ht_prs_preg
├── README.md                 # Overview
├── ht_prs_preg_final.rmd     # R markdown for the analysis
├── ht_prs_preg_final.html    # html generated from rmd file 'ht_prs_sex_final.rmd'
├── articles-functions.R      # Minor R functions for the main analysis
├── select_columns.pl         # Perl script to select columns from tsv files by column name
├── prs_calculations		  # Directory: PRS calculations
	├── prepare_ukb_for_pipeline.R		# Preprocess gwas summaries from ukb
	├── prepare_ega_for_pipeline.R		# Preprocess gwas summaries obtained from collaborators
	├── prs_preecl.json					# Inputs for the FinnGen PRS pipeline
	├── prs_sandbox.wdl					# wdl-file for the FinnGen PRS pipeline
	├── PRS_data_preg3c.csv	  			# Metafile for FinnGen PRS pipeline  
	├── PRS_regions.wdl					# Auxilliary files for FinnGen PRS pipeline (essentially empty file)		

```
