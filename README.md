# dwarf_hamster_hybrids

Hello! This repository contains scripts associated with Hunnicutt et al. 2023. The repository contains the following directories/files of interest:

- **Count_file_generation.md**: This file contains example scripts used to generate the dwarf hamster countfiles used in the manuscript (*e.g.,* Hamsters_psun_14Mar23_singleMapping_featureCounts.csv and associated files).
- **Mouse_hybrid_expression_for_hamster_comparison.Rmd**: This file contains the R code used to generate the mouse analyses included in the script. For the hamster analysis file (Hamster_hybrids_disrupted_expression.Rmd) to run correctly, this must be run **FIRST**!
- **Hamster_hybrids_disrupted_expression.Rmd**: This file contains the R code used to generate the hamster analyses included in the script. For this analysis file to run correctly, the mouse analysis script (Mouse_hybrid_expression_for_hamster_comparison.Rmd) must be run **FIRST**!
- **circos_processing_scripts/**: This directory contains two subfolders (one for mouse and one for hamsters). These folders contain the necessary data and Processing scripts used to generate Supplement Files S7 and S8 (circos style visualization of hypergeometric tests to look for which chromosomes are over- or under-enriched for differentially expressed genes in hybrids). For these scripts to work, you must have Processing installed and have run both the Mouse and Hamster Rmd files.
- **exported_figures/**: This directory is intentionally empty. When the Mouse Rmd file followed by the Hamster Rmd file are run with the parameter "exportFigs" set to TRUE, it will fill with the unedited figures generated for the manuscript.
- **exported_files/**: This directory is intentially empty. When the Mouse Rmd file followed by the Hamster Rmd file are run, files that were generated throughout the course of analysis will be generated (on which other scripts such as those in circos_processing_scripts/ and go_enrichment/ depend).
- **go_enrichment/**: This directory contains the GO enrichment data from DAVID generated in the manuscript. When the Mouse Rmd file followed by the Hamster Rmd file are run, the background gene lists and the DE lists that were used to generated the GO enrichment results will be generated. 
- **induced/**: This directory is intentially empty. When the Mouse Rmd file followed by the Hamster Rmd file are run, files containing the sets of "induced genes" that were generated throughout the course of analysis will be generated.
- **input_data/**: This directory contains various files needed by the Rmd scripts to run including:
	- Hamster_fertility_phenotype_data_7Mar22.csv: raw fertility data of parental species and hybrid dwarf hamsters
	- Hamsters_psun_14Mar23_multiMapping_featureCounts.csv.gz: countfile of dwarf hamster genes; multiMapped reads were allowed for the generation of this file
	- Hamsters_psun_14Mar23_multiMapping_featureCounts.csv.summary: featureCounts metadata file of dwarf hamster gene counts; multiMapped reads were allowed for the generation of this file
	- Hamsters_psun_14Mar23_singleMapping_featureCounts.csv.gz: countfile of dwarf hamster genes; multiMapped reads were NOT allowed for the generation of this file
	- Hamsters_psun_14Mar23_singleMapping_featureCounts.csv.summary: featureCounts metadata file of dwarf hamster gene counts; multiMapped reads were NOT allowed for the generation of this file
	- SingleCell_all_20Jun20_suspenders_featureCounts.csv.gz: countfile of house mouse genes; multiMapped reads were allowed for the generation of this file; this file was also used in the analyses for Hunnicutt et al. (2022) *Evolution*.
	- SingleCell_all_20Jun20_suspenders_featureCounts.csv.summary: featureCounts metadata file of house mouse genes; multiMapped reads were allowed for the generation of this file; this file was also used in the analyses for Hunnicutt et al. (2022) *Evolution*.
	- SingleCell_no-multi_20Jun20_suspenders_featureCounts.csv.gz: countfile of house mouse genes; multiMapped reads were NOT allowed for the generation of this file; this file was also used in the analyses for Hunnicutt et al. (2022) *Evolution*.
	- SingleCell_no-multi_20Jun20_suspenders_featureCounts.csv.summary: featureCounts metadata file of house mouse genes; multiMapped reads were NOT allowed for the generation of this file; this file was also used in the analyses for Hunnicutt et al. (2022) *Evolution*.
	- curated_purity_genes_hamster.csv: this file contains cell type specific genes used for assessing FACS cell population purity in the manuscript
	- ensembl_allgenes_mm38_96_chrs.txt: this file contains the mouse annotations used in the analysis
	- mouse_hamster_orthology_29May23.csv: this file contains the list of hamster/mouse orthologs used in the analysis. This file *can* be re-generated through by running **Hamster_hybrids_disrupted_expression.Rmd**
	- psun_annotations.csv: this file contains the dwarf hamster annotations used in the analysis
	- psun_chrom_sizes.csv: this is a list of sizes of the dwarf hamster chromosomes in the genome used in the analysis 
	- psun_chromosome_scaffold_key.csv: This file is a key connecting the dwarf hamster genome scaffold names with the chromosome numbers they were assigned in this analysis
	- psun_mmus_blast.csv: This file contains dwarf hamster vs. house mouse blast results from Moore et al. (2022) GBE used for orthology assignment
