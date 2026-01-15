peaks=../../Chromatin/Peaks
annot=../../Annotations


#### IAP elements ####

grep IAP ~/Documents/Repeatmasker/mm_rmask_amante.bed > all_IAP.temp
gzcat ${peaks}/musculus_cuttag/H3K4me3_placenta.bed.gz |\
 intersectBed -u -a all_IAP.temp -b - > IAP_mm_K4me3.bed
rm all_IAP.temp


#### L1Md elements ####

grep L1Md ~/Documents/Repeatmasker/mm_rmask_amante.bed > all_L1Md.temp
gzcat ${peaks}/musculus_cuttag/H3K4me3_placenta.bed.gz |\
 intersectBed -u -a all_L1Md.temp -b - > L1Md_mm_K4me3.bed
rm all_L1Md.temp


#### RMER10A elements ####


#get H3K27ac orthologous elements
gzcat ${peaks}/musculus_cuttag/H3K27ac_placenta.bed.gz |\
 intersectBed -u -a ${annot}/Orthology/RMER10A_orth_mm.bed -b - > RMER10A_orth_mm_K27ac.bed
gzcat ${peaks}/pahari_cuttag/H3K27ac_placenta_pahari.bed.gz |\
 intersectBed -u -a ${annot}/Orthology/RMER10A_orth_mp.bed -b - > RMER10A_orth_mp_K27ac.bed

#get epigenetically conserved elements
grep -f <(awk '{print $4}' RMER10A_orth_mp_K27ac.bed) RMER10A_orth_mm_K27ac.bed > RMER10A_epiOrth_mm_K27ac.bed

#get epigenetically divergent elements
grep -vf <(awk '{print $4}' RMER10A_orth_mp_K27ac.bed) RMER10A_orth_mm_K27ac.bed > RMER10A_genOrth_mm_K27ac.bed
grep -vf <(awk '{print $4}' RMER10A_orth_mm_K27ac.bed) RMER10A_orth_mp_K27ac.bed > RMER10A_genOrth_mp_K27ac.bed