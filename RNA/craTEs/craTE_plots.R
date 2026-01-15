library(tidyverse)

load('TS_vs_ES_pseudobulk_tropho_vs_endoth.RData')

volcano = function(data) {
	tib = as_tibble(data$coef, rownames='rep') %>% filter(rep!='(Intercept)')

	ggplot(tib, aes(x=Estimate, y=-log10(p_adj))) +
		geom_point(size=0.5) +
		xlim(-0.75,0.75) +
		ylim(0, 15) +
		xlab('Activity') +
		ylab('-log10 (adj. p-value)') +
		theme_classic()
}

quartz(w=2.5, h=3)
volcano(act_ExE.ectoderm)
volcano(act_SpT)
volcano(act_SynTI)
volcano(act_SynTII)
volcano(act_Glycogen)
volcano(act_S.TGC)

