######################
##APA validation in TCGA
######################
##APA Breast cancer
require(parallel); require(rowr);require(ggplot2)
source('/Users/sinhas8/myCustom_functions.R')
genelist=read.csv('/Users/sinhas8/Downloads/Correlation_and_pvalues.csv')
setwd('/Users/sinhas8/Downloads/')
file_list=list.files()[grep('Combined_PDUIs.txt$',list.files())]
load('/Users/sinhas8/ISLE-Basic/data/TCGA.RData')
prob$race=factor(prob$race)
levels(prob$race)=c('AI', 'AS', 'AA', 'AS', 'HN', 'EA')
more_shortening_in_EA <- function(filename_infunc=file_list[23]){
  cancer_apa=read.csv(filename_infunc, sep='\t')
  # cancer_apa[,1]=sapply(cancer_apa[,1], function(x) strsplit(as.character(x), '\\|')[[1]][2])
  # cancer_apa=cancer_apa[na.omit(match(genelist$gene_name, cancer_apa[,1])),]
  colnames(cancer_apa)=sapply(colnames(cancer_apa), 
                              function(x) paste0(strsplit(x, '\\.')[[1]][1:3], collapse = '-') )
  cancer_apa=cancer_apa[,!is.na(match(colnames(cancer_apa), prob$samples))] 
  dim(cancer_apa)
  infunc_race=prob$race[match(colnames(cancer_apa), prob$samples)]
  infunc_sex=prob$sex[match(colnames(cancer_apa), prob$samples)]
  infunc_age=prob$age[match(colnames(cancer_apa), prob$samples)]
  infunc_stage=prob$stage[match(colnames(cancer_apa), prob$samples)]
  
  cancer_apa=cancer_apa[,which(infunc_race=="AA" | infunc_race=="EA")] 
  infunc_sex=prob$sex[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_age=prob$age[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_stage=prob$stage[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_race=infunc_race[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_race=factor(infunc_race, labels = c('AA', 'EA'))
  return(list(cancer_apa, infunc_race, infunc_sex, infunc_age, infunc_stage, 
              rep(filename_infunc, length(infunc_stage))))
  #Convert PDUi to PSI
  # temp1=cancer_apa
  # cancer_apa= 1-cancer_apa
  # # AA_less=apply(cancer_apa, 1, function(x) err_handle(wilcox.test(x~infunc_race, alternative='l')$p.value) )
  # # c(sum(p.adjust(AA_less, method = 'fdr')<0.05, na.rm = T),
  # #   sum(p.adjust(1-AA_less, method = 'fdr')<0.05, na.rm = T))
  # AA_less=t(apply(cancer_apa, 1, 
  #                 function(x) 
  #                   err_handle(summary(lm(unlist(x)~infunc_race+
  #                                           infunc_age+infunc_sex+infunc_stage))$coefficients['infunc_raceEA',c(1,4)]) ))
  # if(typeof(AA_less)=='list'){
  #   AA_less=do.call(rbind, AA_less)
  # }
  # head(AA_less)
  # c(sum(AA_less[,1]>0 & p.adjust(AA_less[,2], method='fdr')<0.1, na.rm = T),
  #   sum(AA_less[,1]<0 & p.adjust(AA_less[,2], method='fdr')<0.1, na.rm = T))
}

###
apa_dist=mclapply(1:length(file_list),
                  function(x) err_handle(more_shortening_in_EA(file_list[x])), mc.cores = detectCores())
apa_mat_TCGA=do.call(cbind.fill, lapply(apa_dist, function(x) x[[1]] ))
race_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[2]]))) ))
sex_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[3]]))) ))
age_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[4]])) ))
stage_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[5]])) ))
cancertype_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[6]])) ))
sapply(apa_dist[[3]], length)

###
Race=do.call(rbind, sapply(1:length(file_list), 
                           function(x) err_handle(more_shortening_in_EA(file_list[x]))))
df1=data.frame(file_list, APA=do.call(rbind, apa_dist), Race)
df1_corr=data.frame(file_list, APA=do.call(rbind, apa_dist), Race)
apa_dist=mclapply(1:length(file_list), 
                  function(x) err_handle(more_shortening_in_EA(file_list[x])), mc.cores = detectCores())
# Race=do.call(rbind, sapply(1:length(file_list), 
#                            function(x) err_handle(more_shortening_in_EA(file_list[x]))))
df2=data.frame(file_list, APA=do.call(rbind, apa_dist), Race)
df2[df2[,2]>df2[,3],]
df1[df1[,2]>df1[,3],]
df1=na.omit(df1)
write.csv(df1, '/Users/sinhas8/APA_Adriana/TCGA_validation.csv')
######################
##*Checkpoint*
######################
df1[,1]=sapply(df1[,1], function(x) strsplit(as.character(x), '_')[[1]][2] )
colnames(df1)=c('Num of Genes with PSI sig higher in EA', 'Num of Genes with PSI sig higher in AA',
                'AA Samples', 'EA Samples')
######################
##APA diff in different manner
######################
demo=read.csv('/Users/sinhas8/APA_Adriana/2.Data/Demo_Age_Sex.csv')
psi=read.csv('/Users/sinhas8/APA_Adriana/5.Previous_Versions/Finished_Project_APA_prognostic-Marker/Data/PSI_new_imputed.csv')
dim(demo);dim(psi)
colnames(psi)=substring(colnames(psi), 2)
demo=demo[demo$race!='H',]
dim(demo);dim(psi)
psi=psi[,na.omit(match(demo$PATIENT.ACC., colnames(psi)))]
demo$race=factor(as.character(demo$race))

lowPSI_inAA_lm=t(apply(psi, 1, function(x) err_handle( 
  summary(lm(unlist(x) ~ demo$race+demo$gender+demo$age+demo$status))$coefficients['demo$raceEA',c(1,4)] ) ))
head(Gene_PSIhigherinEA)
Gene_PSIhigherinEA=data.frame(GenName=psi[,1][which(lowPSI_inAA_lm[,1]>0 & lowPSI_inAA_lm[,2]<0.1)],
                              lowPSI_inAA_lm[which(lowPSI_inAA_lm[,1]>0 & lowPSI_inAA_lm[,2]<0.1),])
Gene_PSIhigherinAA=data.frame(GenName=psi[,1][which(lowPSI_inAA_lm[,1]<0 & lowPSI_inAA_lm[,2]<0.1)],
                              lowPSI_inAA_lm[which(lowPSI_inAA_lm[,1]<0 & lowPSI_inAA_lm[,2]<0.1),])
write.csv(Gene_PSIhigherinEA, '/Users/sinhas8/APA_Adriana/Gene_PSIhigherinEA_NCIMD.csv')
write.csv(Gene_PSIhigherinAA, '/Users/sinhas8/APA_Adriana/Gene_PSIhigherinAA_NCIMD.csv')

###
psi_all=read.csv('/Users/sinhas8/APA_Adriana/2.Data/PSI.csv')
dim(psi_all)

#########################
##Pan-Cancer TEsting
#########################
more_shortening_in_EA <- function(filename_infunc=file_list[23]){
  cancer_apa=read.csv(filename_infunc, sep='\t')
  # cancer_apa[,1]=sapply(cancer_apa[,1], function(x) strsplit(as.character(x), '\\|')[[1]][2])
  # cancer_apa=cancer_apa[na.omit(match(genelist$gene_name, cancer_apa[,1])),]
  colnames(cancer_apa)=sapply(colnames(cancer_apa), 
                              function(x) paste0(strsplit(x, '\\.')[[1]][1:3], collapse = '-') )
  cancer_apa=cancer_apa[,!is.na(match(colnames(cancer_apa), prob$samples))] 
  col_names=colnames(cancer_apa)
  cancer_apa=t(apply(cancer_apa, 1, scale))
  colnames(cancer_apa)=col_names

  infunc_race=prob$race[match(colnames(cancer_apa), prob$samples)]
  infunc_sex=prob$sex[match(colnames(cancer_apa), prob$samples)]
  infunc_age=prob$age[match(colnames(cancer_apa), prob$samples)]
  infunc_stage=prob$stage[match(colnames(cancer_apa), prob$samples)]
  
  cancer_apa=cancer_apa[,which(infunc_race=="AA" | infunc_race=="EA")] 
  infunc_sex=prob$sex[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_age=prob$age[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_stage=prob$stage[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_race=infunc_race[which(infunc_race=="AA" | infunc_race=="EA")]
  infunc_race=factor(infunc_race, labels = c('AA', 'EA'))
  return(list(cancer_apa, infunc_race, infunc_sex, infunc_age, infunc_stage, 
              rep(filename_infunc, length(infunc_stage))))
  #Convert PDUi to PSI
  # temp1=cancer_apa
  # cancer_apa= 1-cancer_apa
  # # AA_less=apply(cancer_apa, 1, function(x) err_handle(wilcox.test(x~infunc_race, alternative='l')$p.value) )
  # # c(sum(p.adjust(AA_less, method = 'fdr')<0.05, na.rm = T),
  # #   sum(p.adjust(1-AA_less, method = 'fdr')<0.05, na.rm = T))
  # AA_less=t(apply(cancer_apa, 1, 
  #                 function(x) 
  #                   err_handle(summary(lm(unlist(x)~infunc_race+
  #                                           infunc_age+infunc_sex+infunc_stage))$coefficients['infunc_raceEA',c(1,4)]) ))
  # if(typeof(AA_less)=='list'){
  #   AA_less=do.call(rbind, AA_less)
  # }
  # head(AA_less)
  # c(sum(AA_less[,1]>0 & p.adjust(AA_less[,2], method='fdr')<0.1, na.rm = T),
  #   sum(AA_less[,1]<0 & p.adjust(AA_less[,2], method='fdr')<0.1, na.rm = T))
}

###
apa_dist=mclapply(1:length(file_list), 
                  function(x) err_handle(more_shortening_in_EA(file_list[x])), mc.cores = detectCores())
apa_mat_TCGA=do.call(cbind.fill, lapply(apa_dist, function(x) x[[1]] ))
race_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[2]]))) ))
sex_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[3]]))) ))
age_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[4]])) ))
stage_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[5]])) ))
cancertype_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[6]])) ))
psi_tcga= 1-apa_mat_TCGA
psi_tcga
AA_less_lm=do.call(rbind, mclapply(1:nrow(psi_tcga), function(x)
                  err_handle(summary(lm(unlist(psi_tcga[x,])~race_TCGA+
                                          sex_TCGA+age_TCGA+stage_TCGA+cancertype_TCGA))$coefficients['race_TCGAEA',c(1,4)]) ,
                  mc.cores = detectCores()))
c(sum(AA_less_lm[,1]>0 & p.adjust(AA_less_lm[,2], method='fdr')<0.1, na.rm = T),
  sum(AA_less_lm[,1]<0 & p.adjust(AA_less_lm[,2], method='fdr')<0.1, na.rm = T))
FC=unlist(mclapply(1:nrow(psi_tcga), function(x) err_handle( median(unlist(psi_tcga[x,])[race_TCGA=='AA'], na.rm = T)/
                                                               median(unlist(psi_tcga[x,])[race_TCGA=='EA'], na.rm = T) ),
                   mc.cores = detectCores()))

######################
##Plot TCGA results:: Pan_Cancer
######################
FC_thresh=0.13
df2plot=data.frame(Log2Sig= -log(fdrcorr(AA_less_lm[,2]), 10),
                   Log2FC= log(FC, 2),
                   clr=0)
df2plot=df2plot[is.finite(df2plot$Log2FC) & is.finite(df2plot$Log2Sig),]
df2plot$clr[df2plot$Log2Sig>1 & df2plot$Log2FC>  FC_thresh]  =1
df2plot$clr[df2plot$Log2Sig>1 & df2plot$Log2FC< -FC_thresh]=2
write.csv(df2plot, '/Users/sinhas8/APA_Adriana/2.Data/PSI_FC_panCan.csv')
# read.csv(df2plot, '/Users/sinhas8/APA_Adriana/2.Data/PSI_FC_panCan.csv')
FC_thresh=0.13
####
pan_df_ann1=data.frame(label=paste("PSI(EA<AA)=",sum(df2plot$Log2Sig>1 & df2plot$Log2FC< -FC_thresh)), x= -3, y=15)
pan_df_ann2=data.frame(label=paste("PSI(EA>AA)=",sum(df2plot$Log2Sig>1 & df2plot$Log2FC>  FC_thresh)),  x= 3, y=15) 
tiff('/Users/sinhas8/APA_Adriana/pancan_PSI_Nov15.tiff')
ggplot(df2plot, aes(y=Log2Sig, x=Log2FC, color=factor(clr)))+
  geom_point(alpha=0.4, size=1.75)+
  geom_hline(yintercept =c(1), linetype='dashed', color='red')+
  geom_vline(xintercept =c(FC_thresh, -FC_thresh), linetype='dashed', color='blue')+
  labs(x='PSI Log2FoldChange from AA to EA', y='log10FDR')+
  scale_color_manual(values = c('black', 'royalblue', 'Red'))+
  theme_bw(base_size = 15)+
  theme(legend.position="none")+
  geom_label(inherit.aes = FALSE, data=pan_df_ann1, aes(label=label, x= x, y=y))+
  geom_label(inherit.aes = FALSE, data=pan_df_ann2, aes(label=label, x= x, y=y))
dev.off()  
######################
##Plot TCGA results:: Cancer_Type
######################
FC_thresh=0.13
cancer_type_volcano_input<-function(infunc_cancer_type=levels(cancertype_TCGA)[10]){
  infunc_psi_tcga=psi_tcga[,which(cancertype_TCGA==infunc_cancer_type)]
  infunc_race_TCGA=race_TCGA[which(cancertype_TCGA==infunc_cancer_type)]
  infunc_sex_TCGA=sex_TCGA[which(cancertype_TCGA==infunc_cancer_type)]
  infunc_age_TCGA=age_TCGA[which(cancertype_TCGA==infunc_cancer_type)]
  infunc_stage_TCGA=stage_TCGA[which(cancertype_TCGA==infunc_cancer_type)]
  infunc_AA_less_lm=do.call(rbind, mclapply(1:nrow(infunc_psi_tcga), function(x)
    err_handle(summary(lm(unlist(infunc_psi_tcga[x,])~infunc_race_TCGA+
                            infunc_sex_TCGA+infunc_age_TCGA+infunc_stage_TCGA))$coefficients['infunc_race_TCGAEA',c(1,4)]) , mc.cores = detectCores()) )
  infunc_FC=unlist(mclapply(1:nrow(infunc_psi_tcga), function(x) err_handle( median(unlist(infunc_psi_tcga[x,])[infunc_race_TCGA=='AA'], na.rm = T)/
                                                                               median(unlist(infunc_psi_tcga[x,])[infunc_race_TCGA=='EA'], na.rm = T) ),
                            mc.cores=detectCores()))
  return(data.frame(Sig= infunc_AA_less_lm[,2], Log2FC=log(infunc_FC, 2), cancer_type=infunc_cancer_type))
}
length(levels(factor(cancertype_TCGA)))
levels(factor(as.character(df2plot$cancer_type)))
df2plot_cancertype_list=lapply(levels(factor(cancertype_TCGA)), 
                               function(x) err_handle(cancer_type_volcano_input(x)))
df2plot_cancertype=do.call(rbind, df2plot_cancertype_list)
colnames(df2plot_cancertype)[2]='Log2FC'
df2plot_cancertype$clr=0
df2plot_cancertype=df2plot_cancertype[is.finite(df2plot_cancertype$Log2FC) & is.finite(df2plot_cancertype$Sig),]
df2plot_cancertype$Log10Sig = -log(df2plot_cancertype$Sig)
df2plot_cancertype$clr[df2plot_cancertype$Log10Sig > 1 & df2plot_cancertype$Log2FC>  FC_thresh]=1
df2plot_cancertype$clr[df2plot_cancertype$Log10Sig > 1 & df2plot_cancertype$Log2FC< -FC_thresh]=2
# levels(df2plot_cancertype$cancer_type)=sapply(as.character(levels(df2plot_cancertype$cancer_type)), function(x) strsplit(x, '_')[[1]][2])

ann_df1=data.frame(cancer_type=aggregate(df2plot_cancertype$Log2FC, by=list(df2plot_cancertype$cancer_type), function(x) min(x)/2)[,1],
           Log2FC=aggregate(df2plot_cancertype$Log2FC, by=list(df2plot_cancertype$cancer_type), function(x) min(x)/2)[,2],
           Log10Sig=aggregate(df2plot_cancertype$Log10Sig, by=list(df2plot_cancertype$cancer_type), function(x) max(x)/2)[,2],
           Left_score=sapply(split(df2plot_cancertype, df2plot_cancertype$cancer_type), function(x) sum(x[,5]>1 & x[,2]< -FC_thresh ))
           )
ann_df2=data.frame(cancer_type=aggregate(df2plot_cancertype$Log2FC, by=list(df2plot_cancertype$cancer_type), function(x) min(x)/2)[,1],
           Log2FC=aggregate(df2plot_cancertype$Log2FC, by=list(df2plot_cancertype$cancer_type), function(x) max(x)/2)[,2],
           Log10Sig=aggregate(df2plot_cancertype$Log10Sig, by=list(df2plot_cancertype$cancer_type), function(x) max(x)/2)[,2],
           Right_score=sapply(split(df2plot_cancertype, df2plot_cancertype$cancer_type), function(x) sum(x[,5]>1 & x[,2]> FC_thresh ))
)
tiff('/Users/sinhas8/APA_Adriana/test_canType_PSI_NOV15.tiff', height = 1200, width=1200)
ggplot(df2plot_cancertype, aes(y=Log10Sig, x=Log2FC, color=factor(clr)))+
  geom_point(alpha=0.4, size=1.75)+
  geom_hline(yintercept =c(1), linetype='dashed', color='red')+
  geom_vline(xintercept =c(FC_thresh, -FC_thresh), linetype='dashed', color='blue')+
  labs(x='PSI Log2FoldChange from AA to EA', y='log10FDR')+
  scale_color_manual(values = c('black', 'royalblue', 'Red'))+
  theme_bw(base_size = 15)+
  theme(legend.position="none")+
  facet_wrap(cancer_type ~., scales = 'free')+
  geom_label(data=ann_df1, aes(label=Left_score, x = Log2FC, y=Log10Sig), inherit.aes = FALSE)+
  geom_label(data=ann_df2, aes(label=Right_score, x = Log2FC, y=Log10Sig), inherit.aes = FALSE)
dev.off()
