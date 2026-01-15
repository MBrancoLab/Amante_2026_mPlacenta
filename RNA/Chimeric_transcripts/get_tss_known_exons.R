library(tidyverse)
library(GenomicRanges)


##get arguments

args = commandArgs(trailingOnly=TRUE)
gtf.file = args[1]
out.file = args[2]


##read gtf file

gtf = read_tsv(gtf.file, comment='#',
	col_names = c('chr','source','feature','start','end','score','strand','frame','attribute'))
gtf$strand[gtf$strand=='.'] = '*'


##get transcript IDs

split.att = strsplit(gtf$attribute, split='[;\"]')

get.att = function(att) {
	unlist(lapply(split.att, function(x) {
			pos = grep(att, x)
			if (length(pos)>0) {
				return(x[pos + 1])
			} else {
				return('')
			}
		}))
}

gtf = add_column(gtf, transcript = get.att('transcript_id'),
					fpkm = get.att('FPKM'),
					gene = get.att('ref_gene_name'))


##parse exons to GRanges object

gtf.exons = filter(gtf, feature=='exon')

gtf.gr = with(gtf.exons, GRanges(seqnames = chr,
			ranges = IRanges(start=start, end=end),
			strand = strand,
			transcript = transcript))


##get known exons from gencode

ref.exons = read_tsv('gencode.vM25.annotation.gtf', comment='#',
	col_names = c('chr','source','feature','start','end','score','strand','frame','attribute')) %>%
	filter(feature=='exon')

ref.gr = with(ref.exons, GRanges(seqnames = chr,
			ranges = IRanges(start=start, end=end),
			strand = strand))


##filter transcripts overlapping known exons

oe = subsetByOverlaps(gtf.gr, ref.gr, ignore.strand=FALSE)

oe.tx = gtf %>% filter(feature=='transcript', transcript %in% unique(oe$transcript)) %>%
	mutate(tss = case_when(strand=='+' ~ start, strand=='-' ~ end))


##write bed file

tr.bed = oe.tx %>% select(chr, tss, gene, fpkm, strand) %>%
	filter(!is.na(tss)) %>%
	mutate(end=tss+1, .after=tss)
write_tsv(tr.bed, out.file, col_names=FALSE)

