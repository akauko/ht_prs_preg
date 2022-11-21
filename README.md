# ht_prs_preg
Code for the publication: 
Nurkkala, Jouko, et al. "Associations of polygenic risk scores for preeclampsia and blood pressure with hypertensive
disorders of pregnancy." J. Hypertens. (2022) accepted.

* Data: FinnGen https://www.finngen.fi/en
* PRS values were calculuted for FinnGen individuals using PRS-CS pipeline with default settings: https://github.com/getian107/PRScs
* Calculations were performed at the FinnGen PRS pipeline: https://github.com/FINNGEN/CS-PRS-pipeline
* BP GWAS summary from UKBB GWAS v3: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=227859291
* Pre-eclampsia GWAS was based on https://doi.org/10.1038/s41467-020-19733-6, but without FINRISK individuals.

```
ht_prs_preg
├── README.md                 	# Overview
├── ht_prs_preg_final.rmd     	# R markdown for the analysis
├── functions_article.R      	# Minor R functions for the main analysis
├── prs_calculations		# Directory: PRS calculations
	├── prepare_ukb_for_pipeline.R		# Preprocess gwas summaries from ukb
	├── prepare_ega_for_pipeline.R		# Preprocess gwas summaries obtained from collaborators
	├── prs_preecl.json			# Inputs for the FinnGen PRS pipeline
	├── prs_sandbox.wdl			# wdl-file for the FinnGen PRS pipeline (from FinnGen)
	├── PRS_data_preg3c.csv	  		# Metafile for the FinnGen PRS pipeline  
	├── PRS_regions.txt			# Auxilliary file for the FinnGen PRS pipeline (essentially empty file)		

```
