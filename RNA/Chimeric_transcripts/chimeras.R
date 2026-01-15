library(tidyverse)


##read data

tx.files = list.files(path='chimeras')

tx = tibble()
for (f in tx.files) {
	sub = read_tsv(paste('chimeras',f,sep='/'),
					col_names=c('te.chr','te.start','te.end','te.name','te.sw','te.strand',
								'tx.chr','tx.start','tx.end','tx.name','tx.fpkm','tx.strand')) %>%
		mutate(te.id = paste(te.chr,te.start,te.name,sep='_')) %>%
		add_column(sample = gsub('-chimTx.txt','',f))
	tx = rbind(tx,sub)
}


##overview of IAP/L1 transcrips

filter(tx,grepl('IAP',te.name)) %>%
	group_by(sample) %>%
	summarise(n.te=length(unique(te.id)))

filter(tx,grepl('L1Md',te.name)) %>%
	group_by(sample) %>%
	summarise(n.te=length(unique(te.id)))


##select placental chimeras

plac.tes = tx$te.id[grep('Placenta_musculus', tx$sample)]
te.freq = table(plac.tes)
plac.tx = filter(tx, te.id %in% names(te.freq[te.freq>=2])) #TE in at least 2 out of 3 replicates


##placental IAP/L1 transcripts

iap.tx = filter(plac.tx, grepl('IAP',te.name))
l1.tx = filter(plac.tx, grepl('L1Md',te.name), te.end-te.start>5000)


##write output for manual curation

write_tsv(rbind(iap.tx, l1.tx), 'iap_l1_transcripts.txt')

