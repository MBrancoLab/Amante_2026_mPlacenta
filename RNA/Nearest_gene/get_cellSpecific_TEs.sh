peaks=../../Chromatin/Peaks/musculus_cuttag
annot=../../Annotations


#### LTR enhancers ####

#get TSC-specific H3K27ac RLTR13D5
gzcat ${peaks}/H3K27ac_GFPTSCs.bed.gz | intersectBed -u -a ${annot}/RLTR13D5.bed -b - > temp.bed
gzcat ${peaks}/H3K27ac_placenta.bed.gz | intersectBed -v -a temp.bed -b - > TSC_RLTR13D5.bed
rm temp.bed

#get placenta-specific H3K27ac RMER10A
gzcat ${peaks}/H3K27ac_placenta.bed.gz | intersectBed -u -a ${annot}/RMER10A.bed -b - > temp.bed
gzcat ${peaks}/H3K27ac_GFPTSCs.bed.gz | intersectBed -v -a temp.bed -b - > Placenta_RMER10A.bed
rm temp.bed


#### IAP elements ####

grep IAP ~/Documents/Repeatmasker/mm_rmask_amante.bed > all_IAP.temp

#get common H3K4me3 IAPs
gzcat ${peaks}/H3K4me3_GFPTSCs.bed.gz | intersectBed -u -a all_IAP.temp -b - > temp.bed
gzcat ${peaks}/H3K4me3_placenta.bed.gz | intersectBed -u -a temp.bed -b - > Common_IAP.bed
rm temp.bed

#get placenta-specific H3K4me3 IAPs
gzcat ${peaks}/H3K4me3_placenta.bed.gz | intersectBed -u -a all_IAP.temp -b - > temp.bed
gzcat ${peaks}/H3K4me3_GFPTSCs.bed.gz | intersectBed -v -a temp.bed -b - > Placenta_IAP.bed
rm temp.bed all_IAP.temp


#### L1Md elements ####

grep L1Md ~/Documents/Repeatmasker/mm_rmask_amante.bed > all_L1Md.temp

#get common H3K4me3 L1Md
gzcat ${peaks}/H3K4me3_GFPTSCs.bed.gz | intersectBed -u -a all_L1Md.temp -b - > temp.bed
gzcat ${peaks}/H3K4me3_placenta.bed.gz | intersectBed -u -a temp.bed -b - > Common_L1Md.bed
rm temp.bed

#get placenta-specific H3K4me3 L1Md
gzcat ${peaks}/H3K4me3_placenta.bed.gz | intersectBed -u -a all_L1Md.temp -b - > temp.bed
gzcat ${peaks}/H3K4me3_GFPTSCs.bed.gz | intersectBed -v -a temp.bed -b - > Placenta_L1Md.bed
rm temp.bed all_L1Md.temp

