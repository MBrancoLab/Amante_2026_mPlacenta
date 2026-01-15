##Filter repeatmasker annotations and convert to bed files


# musculus mm10, from repeatmasker.org, Dfam 2.0 library
# https://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm406-dfam2.0/mm10.fa.out.gz

# pahari GCF_900095145.1, from UCSC 
# https://hgdownload.soe.ucsc.edu/hubs/GCF/900/095/145/GCF_900095145.1/GCF_900095145.1.repeatMasker.out.gz

# caroli GCF_900094665.2, from UCSC - maybe not used
# https://hgdownload.soe.ucsc.edu/hubs/GCF/900/094/665/GCF_900094665.1/GCF_900094665.1.repeatMasker.out.gz

# rat rn7, from UCSC - maybe not used
# https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.out.gz


##read

mm = read.delim('mm10.fa.out.gz', sep='', skip=3, as.is=T, header=F)
mp = read.delim('GCF_900095145.1.repeatMasker.out.gz', sep='', skip=3, as.is=T, header=F)
mc = read.delim('GCF_900094665.1.repeatMasker.out.gz', sep='', skip=3, as.is=T, header=F)
rn = read.delim('rn7.fa.out.gz', sep='', skip=3, as.is=T, header=F)


##remove simple and low complexity repeats, and tRNAs
##no loose scaffolds

filter.rm = function(rmask) {
	exc = rmask$V11=='Simple_repeat' | rmask$V11=='Low_complexity' | grepl('tRNA',rmask$V10)
	return(rmask[!exc & !grepl('_',rmask$V5),])
}

mm2 = filter.rm(mm)
mp2 = filter.rm(mp)
mc2 = filter.rm(mc)
rn2 = filter.rm(rn)


##make bed file

make.bed = function(rmask) {
	strand = rmask$V9
	strand[strand=='C'] = '-'
	bed = data.frame(chr=rmask$V5, start=rmask$V6, end=rmask$V7,
		name=rmask$V10, score=rmask$V1, strand)
	return(bed)
}

mm.bed = make.bed(mm2)
mp.bed = make.bed(mp2)
mc.bed = make.bed(mc2)
rn.bed = make.bed(rn2)


##write

write.table(mm.bed, 'mm_rmask_amante.bed', sep='\t', quote=F, row.names=F, col.names=F)
write.table(mp.bed, 'mp_rmask_amante.bed', sep='\t', quote=F, row.names=F, col.names=F)
write.table(mc.bed, 'mc_rmask_amante.bed', sep='\t', quote=F, row.names=F, col.names=F)
write.table(rn.bed, 'rn_rmask_amante.bed', sep='\t', quote=F, row.names=F, col.names=F)

