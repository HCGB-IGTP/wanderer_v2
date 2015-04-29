
 __      __                    .___                          
/  \    /  \_____    ____    __| _/___________   ___________ 
\   \/\/   /\__  \  /    \  / __ |/ __ \_  __ \_/ __ \_  __ \
 \        /  / __ \|   |  \/ /_/ \  ___/|  | \/\  ___/|  | \/
  \__/\  /  (____  /___|  /\____ |\___  >__|    \___  >__|   
       \/        \/     \/      \/    \/            \/       




WANDERER DOCUMENTATION

1. AIM
2. CITATION
3. FILE LIST FOR METHYLATION
4. FILE LIST FOR EXPRESSION
5. FREQUENTLY ASKED QUESTIONS
6. VERSION



1. AIM *******************************************************************************************************

Wanderer is a very simple and intuitive web tool allowing real time access and visualization of gene expression and DNA methylation profiles from TCGA data using gene targeted queries. Wanderer is addressed to a broad variety of experimentalists and clinicians without deep bioinformatics skills.\\  

Please feel free to contact us at {adiez,imallona,map@imppc.org}.




2. CITATION **************************************************************************************************

If you find this software useful please consider citing our paper

XXXXX




3. FILE LIST FOR METHYLATION *********************************************************************************

Wanderer generates several images and data files. The naming convention includes the gene, dataset and date and hour of the query.

For instance, the predefined query for METHYLATION produces the following files:


3.1. Wanderer_BRCA1_Mean_methylation_brca_Apr_29_2015_at_170440_CEST.pdf
  A vectorial profile of the average methylation for normals (blue) and tumors (red). The gene location and strand is depicted with an arrow. The CpGs showing statistical differences between normal and tumor are highlighted with an asterisk. Probes within a CpG islands are colored in green.

3.2. Wanderer_BRCA1_Mean_methylation_brca_Apr_29_2015_at_170440_CEST.png
  A raster version of the previous plot (2).

3.3. Wanderer_BRCA1_methylation_brca_Apr_29_2015_at_170440_CEST.pdf
  The profile plot in which the beta values for each sample are linked by lines. The normal's data is in the upper panel (blue) and the tumor's data, in the lower one (red). The gene location and strand is depicted with an arrow. Probes within a CpG islands are colored in green.

3.4. Wanderer_BRCA1_methylation_brca_Apr_29_2015_at_170440_CEST.png
  A raster version of the previous plot (3).

3.5. Wanderer_BRCA1_methylation_brca_Normal_Apr_29_2015_at_170440_CEST.csv
  A comma separated data matrix with the probe name (first column) and a column with the beta values for each of the available normal samples. Sample names are in the first row (header).

3.6. Wanderer_BRCA1_methylation_brca_Tumor_Apr_29_2015_at_170440_CEST.csv
  A comma separated data matrix with the probe name (first column) and a column with the beta values for each of the available tumor samples. Sample names are in the first row (header).

3.7. Wanderer_BRCA1_methylation_brca_annotations_and_statistical_analysis_Apr_29_2015_at_170440_CEST.csv
  A comma separated data matrix with the annotation for each of the probes considered, as well some descriptive analysis.
  The columns correspond to:

  - probe, probe name
  - chr, the chromosome
  - cg_start, the genomic position the CpG starts at
  - cg_end, the genomic position the CpG ends at
  - percentgc, the GC content of the illumina 450k array probe
  - probetype, the type of the illumina 450k array probe 
  - probestart, the genomic position the probe starts at
  - probeend, the genomic position the probe ends at
  - genestart, the genomic position the closest gene to the probe starts at
  - geneend, the genomic position the closest gene to the probe ends at
  - genestrand, the gene strand (negative or positive)
  - ENSEMBL_geneID, the closest gene id at ensembl
  - genebiotype, gene biotype (protein coding, retained intron...)
  - genename, gene symbol
  - cpgistart, for those CpG are inside a CpG island, the coordinate this islands starts at
  - cpgiend, for those CpG are inside a CpG island, the coordinate this islands ends at
  - cpgiid, CpG island identifier
  - Norm_nsamples, number of normal samples for this dataset in this release
  - Norm_mean, mean beta value for normals 
  - Norm_sd, standard deviation for the beta values for normals
  - Tum_nsamples, number of tumor samples for this dataset in this release
  - Tum_mean, mean beta value for tumors
  - Tum_sd, standard deviation for the beta values for tumors
  - wilcox_stat, Wilcoxon Rank Sum Test W parameter (nonparametric comparison of normals vs tumors)
  - pval, Wilcoxon Rank Sum Test p value (nonparametric comparison of normals vs tumors, low values indicates that differences are detected)
  - adj.pval, Benjamini and Hochberg adjustement (False Discovery Rate) for multiple testing.



4. FILE LIST FOR EXPRESSION **********************************************************************************

The predefined query for EXPRESSION produces the following files:


4.1. Wanderer_BRCA1_Mean_expression_brca_Apr_29_2015_at_170440_CEST.pdf
  A vectorial profile of the average methylation for normals (blue) and tumors (red). The gene location and strand is depicted with an arrow. The CpGs showing statistical differences between normal and tumor are highlighted with an asterisk. Probes within a CpG islands are colored in green.


4.2. Wanderer_BRCA1_Mean_expression_brca_Apr_29_2015_at_170440_CEST.png
  A raster version of the previous plot.

Wanderer_BRCA1_boxplot_expression_brca_Apr_29_2015_at_170440_CEST.pdf
Wanderer_BRCA1_boxplot_expression_brca_Apr_29_2015_at_170440_CEST.png
Wanderer_BRCA1_expression_brca_Apr_29_2015_at_170440_CEST.pdf
Wanderer_BRCA1_expression_brca_Apr_29_2015_at_170440_CEST.png
Wanderer_BRCA1_expression_brca_Normal_Apr_29_2015_at_170440_CEST.csv
Wanderer_BRCA1_expression_brca_Normal_RNAseqGENE_Apr_29_2015_at_170440_CEST.csv
Wanderer_BRCA1_expression_brca_Tumor_Apr_29_2015_at_170440_CEST.csv
Wanderer_BRCA1_expression_brca_Tumor_RNAseqGENE_Apr_29_2015_at_170440_CEST.csv
Wanderer_BRCA1_expression_brca_annotations_and_statistical_analysis_Apr_29_2015_at_170440_CEST.csv




5. FREQUENTLY ASKED QUESTIONS  *******************************************************************************

1. Wanderer says there are not enough samples to perform statistical analysis, what does it mean?

2. The CSV files with the data are in a single line!

3. How can I open a CSV in my spreadsheet software?

4. Are you updating this data? When did you take the last snapshot from the TCGA?

5. Which is the annotations origin for gene names and the like?
 Kobor

6. What RPKM means? And RSEM? Why do you use them both?


* VERSION ****************************************************************************************************

Wanderer commit XXX.

This document was generated on April 29th 2015 by Izaskun Mallona.
