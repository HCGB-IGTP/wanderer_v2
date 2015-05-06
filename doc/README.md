# WANDERER DOCUMENTATION ####################################################################################

This is file describes the outputs obtained with TCGA Wanderer.

Please feel free to contact us at adiez imallona or map at imppc.org.

1. AIM
2. CITATION
3. FILE LIST FOR METHYLATION
4. FILE LIST FOR EXPRESSION
5. FREQUENTLY ASKED QUESTIONS
6. VERSION


## AIM #######################################################################################################

Wanderer is a very simple and intuitive web tool allowing real time access and visualization of gene expression and DNA methylation profiles from TCGA data using gene targeted queries. Wanderer is addressed to a broad variety of experimentalists and clinicians without deep bioinformatics skills.


## CITATION ##################################################################################################

If you find this software useful please consider citing our paper (in press).


## FILE LIST FOR METHYLATION #################################################################################

Wanderer generates several images and data files. The naming convention includes the gene, dataset and date and hour of the query.

For instance, the predefined query for METHYLATION produces the following files:

1. `Wanderer_BRCA1_methylation_brca_Apr_29_2015_at_170440_CEST.pdf`

  A vectorial profile plot in which the beta values for each sample are linked by lines. The normal's data is in the upper panel (blue) and the tumor's data, in the lower one (red). The gene location and strand is depicted with an arrow. Probes labels within a CpG islands are colored in green.

2. `Wanderer_BRCA1_methylation_brca_Apr_29_2015_at_170440_CEST.png`

  A raster version of the previous plot.

3. `Wanderer_BRCA1_Mean_methylation_brca_Apr_29_2015_at_170440_CEST.pdf`

  A vectorial profile plot of the average methylation for normals (blue) and tumors (red). The gene location and strand is depicted with an arrow. The CpGs showing statistical differences between normal and tumor are highlighted with an asterisk (Wilcoxon adjusted p-value < 0.05). Probes labels within a CpG islands are colored in green.

4. `Wanderer_BRCA1_Mean_methylation_brca_Apr_29_2015_at_170440_CEST.png`

  A raster version of the previous plot.

5. `Wanderer_BRCA1_methylation_brca_Normal_Apr_29_2015_at_170440_CEST.csv`

  A comma separated data matrix with the probe name (first column) and a column with the beta values for each of the available normal samples. Sample names are in the first row (header).

6. `Wanderer_BRCA1_methylation_brca_Tumor_Apr_29_2015_at_170440_CEST.csv`

  A comma separated data matrix with the probe name (first column) and a column with the beta values for each of the available tumor samples. Sample names are in the first row (header).

7. `Wanderer_BRCA1_methylation_brca_annotations_and_statistical_analysis_Apr_29_2015_at_170440_CEST.csv`

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
  - genestrand, the closest gene strand
  - ENSEMBL_geneID, the closest gene id at ensembl
  - genebiotype, the closest gene biotype (protein coding, retained intron...)
  - genename, the closest gene symbol
  - cpgistart, for those CpG are inside a CpG island, the coordinate this islands starts at
  - cpgiend, for those CpG are inside a CpG island, the coordinate this islands ends at
  - cpgiid, for those CpG are inside a CpG island, CpG island identifier
  - Norm_nsamples, number of normal samples for this dataset in this release
  - Norm_mean, mean of the beta values for normals 
  - Norm_sd, standard deviation of the beta values for normals
  - Tum_nsamples, number of tumor samples for this dataset in this release
  - Tum_mean, mean of the beta values for tumors
  - Tum_sd, standard deviation of the beta values for tumors
  - wilcox_stat, Wilcoxon Rank Sum Test W parameter (nonparametric comparison of normals vs tumors)
  - pval, Wilcoxon Rank Sum Test p value (nonparametric comparison of normals vs tumors, low values indicates that differences are detected)
  - adj.pval, Benjamini and Hochberg adjustement (False Discovery Rate) for multiple testing.



## FILE LIST FOR EXPRESSION ##################################################################################

The predefined query for EXPRESSION produces the following files:

1. `Wanderer_BRCA1_expression_brca_Apr_29_2015_at_170440_CEST.pdf`

  A vectorial profile plot in which the log2-transformed RPKM values for each exon and each sample are linked by lines. The normal's data is in the upper panel (blue) and the tumor's data, in the lower one (red). The gene location and strand is depicted with an arrow.

2. `Wanderer_BRCA1_expression_brca_Apr_29_2015_at_170440_CEST.png`

  A raster version of the previous plot.

3. `Wanderer_BRCA1_Mean_expression_brca_Apr_29_2015_at_170440_CEST.pdf`

  A vectorial profile plot of the average expression for normals (blue) and tumors (red). The gene location and strand is depicted with an arrow. The exons showing statistical differences between normal and tumor are highlighted with an asterisk (Wilcoxon adjusted p-value < 0.05).

4. `Wanderer_BRCA1_Mean_expression_brca_Apr_29_2015_at_170440_CEST.png`

  A raster version of the previous plot.

5. Wanderer_BRCA1_boxplot_expression_brca_Apr_29_2015_at_170440_CEST.pdf`

  A vectorial boxplot and a stripchart showing the expression values summarized by gene for normals (blue) and tumors (red). The expression values are log2-transformed normalized RSEM values and reflect the expression of the gene as a whole.

6. `Wanderer_BRCA1_boxplot_expression_brca_Apr_29_2015_at_170440_CEST.png`

  A raster version of the previous plot.
 
7. `Wanderer_BRCA1_expression_brca_Normal_Apr_29_2015_at_170440_CEST.csv`

  A comma separated data matrix with the exon name (first column) and a column with the log2-transformed RPKM values for each of the exons of every available normal sample. Sample names are in the first row (header).

8. `Wanderer_BRCA1_expression_brca_Normal_RNAseqGENE_Apr_29_2015_at_170440_CEST.csv`

   A comma separated data matrix with the available normal sample names (first column) and a column with the log2-transformed RSEM values for your gene of interest.

9. `Wanderer_BRCA1_expression_brca_Tumor_Apr_29_2015_at_170440_CEST.csv`

   A comma separated data matrix with the exon name (first column) and a column with the log2-transformed RPKM values for each of the exons of every available tumor sample. Sample names are in the first row (header).

10. `Wanderer_BRCA1_expression_brca_Tumor_RNAseqGENE_Apr_29_2015_at_170440_CEST.csv`

   A comma separated data matrix with the available tumor sample names (first column) and a column with the log2-transformed RSEM value for your gene of interest.

11. `Wanderer_BRCA1_expression_brca_annotations_and_statistical_analysis_Apr_29_2015_at_170440_CEST.csv`

  A comma separated data matrix with the annotation of your gene expression data, as well some descriptive analysis.
  The columns correspond to:

  - exon, exon identifier according to the TCGA pipeline
  - id, exon identifier according to Genome Browser exons track
  - ENSEMBL_geneID, the gene identifier according to ENSEMBL (ENSG identifier)
  - ENSEMBL_transcriptID, the transcript identifier according to ENSEMBL (ENST identifier)
  - chr, the chromosome the exon is located at
  - exon_start, the genomic coordinate the exon starts at 
  - exon_end, the genomic coordinate the exon ends at
  - strand, the genomic strand of the exon's gene.
  - genestart, the genomic start position of the exon's gene.
  - geneend, the genomic end position of the exon's gene.
  - genebiotype, exon's gene biotype (protein coding, retained intron...)
  - genename, exon's gene symbol
  - rnaseqgeneid, the TCGA exon's gene identifier
  - Norm_nsamples, number of normal samples for this dataset in this release
  - Norm_mean, mean of log2-transformed RPKM for normals 
  - Norm_sd, standard deviation of log2-transformed RPKM in normals
  - Tum_nsamples, number of tumor samples for this dataset in this release
  - Tum_mean, mean of log2-transformed RPKM for tumors
  - Tum_sd, standard deviation of log2-transformed RPKM in tumors
  - wilcox_stat, Wilcoxon Rank Sum Test W parameter (nonparametric comparison of normals vs tumors)
  - pval, Wilcoxon Rank Sum Test p value (nonparametric comparison of normals vs tumors, low values indicates that differences are detected)
  - adj.pval, Benjamini and Hochberg adjustement (False Discovery Rate) for multiple testing.



## FREQUENTLY ASKED QUESTIONS  ###############################################################################

1. Wanderer says there are not enough samples to perform statistical analysis, what does it mean?

  We compute the Wilcoxon test given that there are at least a normal and a tumor samples; if no tumors or no normals are available, we cannot test whether they are equivalent or not.
  We note that, although we calculate the test with small number of samples, this result might be meaningless without a minimum number of cases.

2. The CSV files are a single line file!

  We use Linux/UNIX linefeeds. Although your text viewer (i.e. notepad) might merge all the lines under Windows, your spreadsheet software (i.e. Excel) will recognise the line separation properly.

3. How can I open a CSV in my spreadsheet software?

  Just use the 'open with' (right click with the mouse). Or just lauch the application and import data from text, comma-separated or character-delimited values (this depends on the software).

4. How and when did you download the data from the TCGA?

   Data from TCGA was dowloaded using the TCGA-Assembler with the Traverse Result File 'DirectoryTraverseResult_Jul-08-2014.rda'
Y. Zhu, P. Qiu, Y. Ji*. TCGA-Assembler: Open-Source Software for Retrieving and Processing TCGA Data. Nature Methods. 11:599-600, 2014. | doi:10.1038/nmeth.2956)

5. Which is the annotations origin for gene names and the like?

  For 450k methylation array data we use the nice paper from Price, Kobor et al. named _Additional annotation enhances potential for biologically-relevant analysis of the Illumina Infinium HumanMethylation450 BeadChip array_ Epigenetics & Chromatin 2013, 6:4  doi:10.1186/1756-8935-6-4.

  For RNAseq data, since TCGA only offers the position of the exons but not the gene it is associated to, we find the associated gene and some other annotation information merging BioMart and Ensembl gene annotations with the exon position comming from TCGA and Genome Browser. To merge this data we use the function overlapRegions from the Bioconductor package called regioneR.


6. What does RPKM mean? And RSEM? Why do you use them?

  We take advantage of the TCGA pipeline for RNA seq v2 expression analysis. RPKM are used for exons and RSEM for genes. We log2-transform this values adding one to get rid of zeroes; that means we plot `log_2(x+1)` RPKM or RSEM values.



## VERSION ####################################################################################################

`
commit 51c3d0597cc7e83521bff236f5a26d187e71870b
Author: Izaskun Mallona <imallona@imppc.org>
Date:   Tue Apr 28 12:20:53 2015 +0200
`