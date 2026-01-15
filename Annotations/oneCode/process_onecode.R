

##read onecode output; ignore headers

raw = scan('All_IAPs.txt.gz', character(), comment='S', sep='\n')


##parse grouped annotations into list

ltr.head = grep('#', raw)

ltr = list()
for (i in 1:length(ltr.head)) {
	start = ltr.head[i] + 1
	if (i==length(ltr.head)) {
		end = length(ltr.head)
	} else {
		end = ltr.head[i+1] - 1
	}
	
	if (end < start) {
		#single entries
		entry = strsplit(gsub('#','',raw[start-1]), split='\t')
	} else {
		#grouped entries
		entry = strsplit(raw[start:end], split='\t')
	}
	
	te = unlist(lapply(entry, function(x) x[10]))
	chr = unlist(lapply(entry, function(x) x[5]))
	start = unlist(lapply(entry, function(x) as.numeric(x[6])))
	end = unlist(lapply(entry, function(x) as.numeric(x[7])))
	length = unlist(lapply(entry, function(x) as.numeric(x[8])))
	strand = unlist(lapply(entry, function(x) x[9]))
	
	ltr[[i]] = data.frame(te, chr, start, end, length, strand)
}


##classify groups

ltr.class = lapply(ltr, function(x) {
	int.pos = grepl('int',x$te) | grepl('_I',x$te)
	ltr.pos = grepl('LTR',x$te) & !int.pos
	if (any(ltr.pos)) {
		if (all(ltr.pos)) {
			return('solo_LTR')
		} else if (ltr.pos[1] & ltr.pos[nrow(x)]) {
			return('proviral')
		} else {
			return('truncated')
		}
	} else {
		return('int_fragment')
	}
})


##make dataframe

df = data.frame(chr = unlist(lapply(ltr, function(x) x$chr[1])),
			start = unlist(lapply(ltr, function(x) x$start[1])),
			end = unlist(lapply(ltr, function(x) x$end[nrow(x)])),
			strand = unlist(lapply(ltr, function(x) x$strand[1])),
			type = unlist(ltr.class),
			repname = unlist(lapply(ltr, function(x) x$te[1])),
			structure = unlist(lapply(ltr, function(x) paste(x$te, collapse=':'))))
df$length = df$end - df$start
df$strand[df$strand=='C'] = '-'


##write all merged elements

write.table(df, 'IAP_merged_mm10.txt', sep='\t', quote=F, row.names=F)


##write bed file of proviral elements

write.table(df[df$type=='proviral',c(1:3,6,8,4)], 'IAP_proviral_mm10.bed', sep='\t', quote=F, row.names=F, col.names=F)

