# PSI density for breast and Pan-Cancer.
require(ggplot2)
apa_dist=mclapply(1:length(file_list), 
                  function(x) err_handle(more_shortening_in_EA(file_list[x])), mc.cores = detectCores())
apa_mat_TCGA=do.call(cbind.fill, lapply(apa_dist, function(x) x[[1]] ))
race_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[2]]))) ))
sex_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[3]]))) ))
age_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[4]])) ))
stage_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[5]])) ))
cancertype_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[6]])) ))
cancertype_TCGA=factor(cancertype_TCGA)
levels(cancertype_TCGA)=sapply(as.character(levels(cancertype_TCGA)), function(x) strsplit(x, '_')[[1]][2])
##############################################
#PanCan
##############################################
dat=data.frame(PSI=c(apply(apa_mat_TCGA[,which(cancertype_TCGA=='BRCA')], 1, function(x) x)),
               Hist=cancertype_TCGA[which(cancertype_TCGA=='BRCA')],
               Race=race_TCGA[which(cancertype_TCGA=='BRCA')])
tiff('/Users/sinhas8/APA_Adriana/5.Previous_Versions/Finished_Project_APA_prognostic-Marker/Results/BRCA_Tumor_PSIDensity_byRace.tiff', 
     width=600, height=600)
ggplot(dat, aes(x=PSI, fill=Race, colour=Race))+
  geom_density(alpha=0.3, position="identity", linetype='dashed', size=0.1)+
  theme_bw(base_size = 20)+
  scale_fill_manual(values = c('#e41a1c','#377eb8'))
dev.off()

dat=data.frame(PSI=c(apply(apa_mat_TCGA, 1, function(x) x)),
               Hist=cancertype_TCGA,
               Race=race_TCGA)
rm(prob)
tiff('/Users/sinhas8/APA_Adriana/5.Previous_Versions/Finished_Project_APA_prognostic-Marker/Results/PanCan_Tumor_PSIDensity_byRace.tiff', 
     width=600, height=600)
ggplot(dat, aes(x=PSI, fill=Race, colour=Race))+
  geom_density(alpha=0.3, position="identity", linetype='dashed', size=0.1)+
  theme_bw(base_size = 20)+
  scale_fill_manual(values = c('#e41a1c','#377eb8'))
dev.off()

# For all cancer type
tiff('/Users/sinhas8/APA_Adriana/5.Previous_Versions/Finished_Project_APA_prognostic-Marker/Results/For_allCancerTypes_Tumor_PSIDensity_byRace.tiff', 
     width=2400, height=2400)
ggplot(dat, aes(x=PSI, fill=Race, colour=Race))+
  geom_density(alpha=0.3, position="identity", linetype='dashed', size=0.1)+
  facet_wrap(.~Hist)+
  theme_bw(base_size = 25)+
  scale_fill_manual(values = c('#e41a1c','#377eb8'))
dev.off()

# Write Table 1
write.csv(Demo_Age_sex[,-c(1, 5)], '/Users/sinhas8/APA_Adriana/Resubmission Main Files/Revised_Table_1.csv')