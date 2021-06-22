# Identify whether PSI load affects survival
require(survival)
require(survminer)
median(colMedians(as.matrix(PSI_Norm)))
delta_PSI_load=colMedians(as.matrix(delta_PSI))
dat=data.frame(delta_PSI=delta_PSI_load,
               Hist=histology_df,
               Race=race_df,
               demo[match(gsub('X','',names(delta_PSI_load)), demo$acc),]
               )
dat$status=as.numeric(as.character(factor(dat$status, labels = c('0', '1'))))
dat$delta_PSI_quantized= xtile(dat$delta_PSI, 4)
dat=dat[dat$delta_PSI_quantized==1 | dat$delta_PSI_quantized==4,]
dat=na.omit(dat)
fit=survfit(Surv(survival,  dat$) ~ delta_PSI_quantized, data=dat )
ggsurvplot(fit, data=dat, pval = TRUE, pval.method = TRUE)
