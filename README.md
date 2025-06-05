# M1-Internship
M1 internship, april-june 2025, in Habermann Lab (Computational Biology)


Project title: Friedreich's Ataxia, RNA-Seq data analysis and metabolic modelling

Project description: 
1. Data
Bulk-RNA sequencing in cardiac cells of WT and KO mice at three weeks (WT3w and KO3w, respectively), and WT, KO-NaCl and KO treated with AAV mice at eight weeks (WT8w-NaCl, KO8w-NaCl and KO8w-AAV, respectively).

2. Differential gene expression analyses
With DESeq2 package (version 1.46.0)

3. Gene selection for ataxiaXplorer
FRDA gene list built containing genes that were identified from differential expression analyses of FRDA data, that were not found in the ataxia interactome builyt by the lab from SCA ataxia data.

4. Functional enrichment analyses and annotation
Enrichment analysis, using Enrichr, was performed on the DEGs identified for each condition comparison: WT8w-NaCl versus KO8w-NaCl, KO8w-NaCl versus KO8w-AAV, KO8w-AAV versus WT8w-NaCl and WT3w versus KO3w.

5. Metabolic modelling
Metabolic modelling of the WT and KO phenotypes at three weeks was performed using mitoMAMMAL (pFBA).
