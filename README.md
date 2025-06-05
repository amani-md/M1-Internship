# M1-Internship

M1 internship, April–June 2025, in the Habermann Lab (Computational Biology)

## Project Title
Friedreich's Ataxia: RNA-Seq Data Analysis and Metabolic Modelling

## Project Description

1. **Data**  
   Bulk RNA sequencing was performed on cardiac cells from wild-type (WT) and knockout (KO) mice at three weeks (WT3w and KO3w, respectively), and from WT, KO-NaCl, and KO treated with AAV at eight weeks (WT8w-NaCl, KO8w-NaCl, and KO8w-AAV).

2. **Differential Gene Expression Analyses**  
   Conducted using the DESeq2 package (version 1.46.0).

3. **Gene Selection for ataxiaXplorer**  
   An FRDA gene list was created, including genes identified from differential expression analyses of FRDA data that were not present in the ataxia interactome previously built by the lab using SCA ataxia data.

4. **Functional Enrichment Analyses and Annotation**  
   Enrichment analysis (using Enrichr) was performed on the DEGs identified in the following comparisons:  
   - WT8w-NaCl vs KO8w-NaCl  
   - KO8w-NaCl vs KO8w-AAV  
   - KO8w-AAV vs WT8w-NaCl  
   - WT3w vs KO3w

5. **Metabolic Modelling**  
   Performed using mitoMAMMAL (pFBA) to model the metabolism of WT and KO mice at three weeks.

## Repository Layout

1. **Statistical_Analysis/**  
   Contains all R Markdown files used for exploratory analysis, differential expression, enrichment, and annotation. Also includes the resulting figures and data outputs (e.g., DEG lists, FRDA gene list, and DESeq2 results for all genes in all comparisons).

2. **Metabolic_Modeling/**  
   Contains the Jupyter notebook used for metabolic modeling, associated data files (e.g., modified read counts, mitoMAMMAL files), and resulting figures and output files (predicted fluxes).
