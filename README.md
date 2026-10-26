# Amante_2026_mPlacenta
Scripts and files used in Amante et al (2026) "IAP retrotransposons contribute to the transcriptional diversity of the murine placenta"

## Primary data processing

*Chromatin (CUT&Tag, ATAC-seq)*

Reads were mapped with bowtie2 (bowtie_PE.sh) and blacklisted regions (generated using get_blacklist.py) filtered out (blacklist_filter.sh). Bam files were converted to bigwig (bam2bigwig.sh). CUT&Tag peaks were detected using GoPeaks (gopeaks.sh) and ATAC-seq peaks with MACS2 (atac-macs2.sh).

*RNA*

Bulk RNA-seq was mapped with HiSat2 and multi-mapping reads filtered out (hisat_PE.sh). Read counts were extracted using Seqmonk. Bam files were converted to bigwig (bam2bigwig.sh). For the analysis of TE expression, reads were mapped with STAR and read counts extracted using SQuIRE and TEtranscripts (squire_tetranscripts.sh). Published scRNA-seq data (downloaded with fastq_dump.sh) was processed with CellRanger (cellranger.sh), generating the matrix.mtx, features.tsv and barcodes.tsv files. A Seurat object was generated from these data, where data were filtered and clusters identified (Marsh_Seurat.R).

## Annotations

Repeatmasker outputs for each species were filtered and parsed to bed format (parse_repeatmasker.R). The 'One Code To Find Them All' tool was used to merge LTR and internal ERV annotations (onecode.sh), with further processing using the 'process_onecode.R' script. For identification of orthologous elements see 'Comparative analyses' section below.

## Cell identity

RT-qPCR data from TSC differentiation were plotted (Figure S1A) using the 'RT-qPCR.R' script. For cell-type deconvolution, Cybersort was used with placental scRNA-seq data and the output processed (cybersort.R) to generate Figure S1B. 

## TE-derived enhancers

*Peak enrichment*

A permutation analysis was performed to identify TE subfamilies enriched for a given set of CUT&Tag or ATAC-seq peaks. Species-specific input files (e.g., mm10_150PE_lowMap.bed, mm10.chrom.sizes) and arguments list (e.g., mm_permutation_args.txt) were fed into the permutation pipeline (permutation_pipeline.R, run via permutation.sh), yielding the files ending in '_overlap.txt' that are grouped in different folders. These were used to generate: 1) Figures 1B,C and S1C with the 'musculus.R' script, 2) Figures S1F-H with the 'encode.R' script. The overall comparison between ENCODE and own data in Figure S1E was generated with the 'source_comparison.R' script in the Peaks folder.

*craTEs*

Pseudobulk counts from trophoblast cell types (in scRNA_processing folder) and early epiblast/extraembryonic ectoderm (make_pseudobulk.R) were used to run craTEs. It's output (TS_vs_ES_pseudobulk_tropho_vs_endoth.RData) was used by the craTE_plots.R script to generate Figure S1D.

*Profiles*

H3K27ac CUT&Tag profiles over selected TE subfamilies (Figure 1D) were generated with deepTools (heatmaps.sh).

*Nearest gene*

Relative gene expression (Placenta over TSC) was coupled to the nearest TE from selected sets using the 'proximity_analysis.R' script. This generated Figure 1E.

## IAP expression

*Profiles*

Mean H3K4me3 signal profiles were generated with deepTools (trends.sh), yielding the outputs in the Trends folder, which were used in the 'trends.R' script to generate Figures 2A and S2A-C.

*Subfamily expression*

Outputs from SQuIRE and TEtranscripts were used in the 'subfamily.R' script to generate Figures 2B and S2D.

*Element expression*

RNA profiles over selected TE subfamilies (Figure S2E) were generated with deepTools (rna_heatmaps.sh). RPM values from uniquely mapped reads were used in the 'IAP_expression.R' scripts to generate Figures 2C,D and S2A. IAP gtf annotations were generated (make_IAP_gtf.R) to append to gene annotations used for scRNA-seq data processing. The respective Seurat object was used in the 'IAP_scRNAseq.R' script to generate Figures 2F,G.

*DNA methylation*

Published placenta and brain DNA methylation data were used in the IAPEz_LTR_WGBS.R to generate Figure S3D. Pyrosequencing data from matched placenta and tail samples were used in the VM-IAP_pyro.R script to generate Figure S3E.

## Chimeric transcripts

Correlations between IAPs and gene expression (Figures 3A,B and S4A,B) were tested with the 'proximity_analysis.R' script in the Nearest_gene folder. Chimeric transcripts were identified by running Stringtie, selecting transcripts overlapping known exons, and intersecting TSSs with the Repeatmasker annotation (stringtie.sh, which also uses get_tss_known_exons.R). Further filtering was performed (chimeras.R) before manual curation. Isoform-specific RT-qPCR data were plotted (Figure 3D) using the 'RT-qPCR.R' script. Transcript expression values from Stringtie (merged_quant folder) were used in the merged_gtf_quant.R' script to generate Figure 3E. The scRNA-seq Seurat object was used in the 'chimeras_scRNAseq.R' script to generate Figure 3F.

## CRISPRi

*sgRNA hits*

The sgRNA sequence was mapped to mm10 using bowtie2 (-a option) and intersected with the Repeatmasker annotation using bedtools (-wao option). The output was used in the 'sgRNA_IAP_hits.R' script to generate Figure 4B.

*Expression changes*

SQuIRE outputs were used in the 'subfamily_expression.R' script to plot Figure S6A. To test effects on gene expression, raw read counts were used: 1) in the 'proximity_analysis.R' script to plot expression changes for different subgroups of proximal IAPs (Figures 4D and S6B), 2) in the 'IAP_CRISPRi_DEseq.R' script to perform differential expression (Figure S6C) and gene ontology (Figure S6D) analyses.

## Comparative analyses

*Orthologous regions*

Orthologous TEs between musculus and other species were identified using LiftOver, then counted and plotted in Figure 5A (liftover.sh and orthology.R scripts in Annotations/Orthology). LiftOver was also used to identify orthologous and non-orthologous CUT&Tag peaks between musculus and pahari, and bedtools used to further subdivide orthologous regions into epigenetically conserved or divergent (peak_orthology.sh, in Peaks folder).  Read counts over orthologous RMER10A elements were processed with the 'RMER10A_orth_counts.R' script to generate Figure 5B.

*Peak enrichment*

Orthologous (epigenetically conserved or divergent) and non-orthologous CUT&Tag peaks were analysed using a permutation pipeline as above (TE-derived enhancers section) to identify enriched TE subfamilies for each group. Outputs from the permutation pipeline were processed using: 1) the 'pahari.R' script to generate Figure S8C, 2) the 'orthology.R' script to generate Figures 5D and S8D,E.

*Gene expression*

Differential expression analysis between musculus and pahari was performed and the plot in Figure S8A generated with the 'placenta_DEseq.R' script. The comparison between musculus and pahari tissues (Figure S8B,G) was performed with the 'pahari_vs_encode.R' script. Correlations between species-specific CUT&Tag peaks and gene expression were analysed using the 'proximity_analysis.R' script (Pahari folder) to generate Figures 5E and S8F. For Figure 5E, lists of H3K4me3-marked TEs were generated by the 'get_species_TEs.sh' script.

*Inter-strain polymorphisms*

Lists of structural variants across mouse strains involving IAP elements were used by the 'Ferraj_SVs.R' script (Annotations folder) to generate Figure S9.
