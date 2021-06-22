# APA machinary across Race in TCGA
######################
# Load and Preprocess files
######################
source('/Users/sinhas8/myCustom_functions.R')
require(parallel)
require(ggpubr)
require(ggplot2)
require(matrixStats)
prob_wdMut=readRDS('/Users/sinhas8/Downloads/TCGA_withMut.RDS')
load('/Users/sinhas8/ISLE-Basic/data/TCGA.RData')
APA_mach=read.csv('/Users/sinhas8/Downloads/APA_geneset.csv')
APA_mach_id=match(unlist(APA_mach$GENE.NAME), prob$genes)
CORE_APA_mach=read.csv('/Users/sinhas8/Downloads/geneset (3).txt')
CORE_APA_mach_id=na.omit(match(CORE_APA_mach$Polyadenylation, prob$genes))
prob$race=factor(prob$race)
levels(prob$race)=c('AI', 'AS', 'AA', 'AS', 'HN', 'EA')
test_machinary_activity_acrossRace<-function(infunc_cancer_type=levels(factor(prob$types))[1], 
                                             Score=mach_act,
                                             hypo='g'){
  wilcox.test(Score[prob$types==infunc_cancer_type & prob$race=='AA'], Score[prob$types==infunc_cancer_type & prob$race=='EA'], alternative = hypo)$p.value
}
mach_act_Exp=colMedians(prob$mRNA[APA_mach_id,])
mach_act_CNV=colMedians(prob$scna[APA_mach_id,])
#mach_act_CNV=colMedians(prob$[APA_mach_id,])
dfsig=data.frame(sig_diff_exp=sapply(levels(factor(prob$types)),function(x) test_machinary_activity_acrossRace(infunc_cancer_type = x, Score = mach_act_Exp, hypo='l')),
           sig_diff_cnv=sapply(levels(factor(prob$types)),function(x) test_machinary_activity_acrossRace(infunc_cancer_type = x, Score = mach_act_CNV, hypo='l'))
)
write.csv(dfsig, '/Users/sinhas8/APA_Adriana/Significance_ofAPAMachinarywithhighActivityinEA.csv')
# cor.test(ann_df1$Left_score -ann_df2$Right_score, ann_df1$Left_score )
######################
##APA machinary activity with PSI
######################
df2plot=read.csv('/Users/sinhas8/Downloads/2.Data/TCGA_PSI.csv')
mach_act_Exp=colMedians(prob_wdMut$mRNA[APA_mach_id,])
df2plot=na.omit(df2plot)
df2plot$APA_machinary_activity=mach_act_Exp[match(df2plot$sample_id, gsub('-','.',colnames(prob_wdMut$Mut)))]
# cor.test(df2plot$PSI, df2plot$APA_machinary_activity, method='spear')
colorType_Set='Set1'
write.csv(df2plot, '/Users/sinhas8/APA_Adriana/2.Data/SuppTable17.csv')
tiff('/Users/sinhas8/APA_Adriana/APAmachinary_byCancerType_Feb6.tif',width = 1800, height = 650)
ggplot(df2plot, aes(x=cancer_type,y=APA_machinary_activity, fill=race))+
  geom_boxplot(data=df2plot, aes(fill=race))+
  labs(title="APA_mach_activity in 23 cancer types (TCGA)",x="Race", y = "APA machinary activity")+
  facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 20)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4)+
  scale_fill_brewer(palette=colorType_Set)
dev.off()
######################
##Corr of CORE APA machinary activity with PSI
######################
mach_act_Exp=colMeans(prob_wdMut$mRNA[CORE_APA_mach_id,])
df2plot$APA_machinary_activity=mach_act_Exp[match(df2plot$sample_id, gsub('-','.',colnames(prob_wdMut$Mut)))]
df2plot=na.omit(df2plot)
#**#
cor.test(df2plot$PSI, df2plot$APA_machinary_activity, method='spear')
tiff('/Users/sinhas8/APA_Adriana/CORE_APAmachinary_byCancerType_Feb6.tif',width = 1800, height = 650)
ggplot(df2plot, aes(x=cancer_type,y=APA_machinary_activity, fill=race))+
  geom_boxplot(data=df2plot, aes(fill=race))+
  labs(title="APA_mach_activity in 23 cancer types (TCGA)",x="Race", y = "APA_machinary_activity")+
  facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 20)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)
dev.off()
######################
##Corr of PSI with mutation states of each gene
######################
df2plot_FF=df2plot[!is.na(match(df2plot$sample_id, gsub('-','.',colnames(prob_wdMut$Mut)))),]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}
df2plot_FF$Normalized_PSI=scaling_cancerType(df2plot_FF$PSI, df2plot_FF$cancer_type)

Mut_forPSI=prob_wdMut$Mut[,match(df2plot_FF$sample_id, gsub('-','.',colnames(prob_wdMut$Mut)))]

cosmic=read.csv('/Users/sinhas8/APA_Adriana/COSMIC.tsv', sep='\t')
drivers_id=na.omit(match(cosmic$Gene.Symbol, rownames(Mut_forPSI)))

driver_effect_onPSI=mclapply(drivers_id, function(x) err_handle(c(p_value=wilcox.test(df2plot_FF$Normalized_PSI ~ unlist(Mut_forPSI[x,]))$p.value,
                                                              FC=median(df2plot_FF$Normalized_PSI[as.logical(Mut_forPSI[x,])], na.rm=T)/
                                                                median(df2plot_FF$Normalized_PSI[!as.logical(Mut_forPSI[x,])], na.rm=T),
                                                              corr_pvalue=summary(lm(df2plot_FF$Normalized_PSI ~ unlist(Mut_forPSI[x,])+df2plot_FF$cancer_type))$coefficients[2,4] )),
                             mc.cores = detectCores())
######################
# ANALYZE RESULTS: FOR ALL DRIVER GENES
######################
names(driver_effect_onPSI)=rownames(Mut_forPSI)[na.omit(match(cosmic$Gene.Symbol, rownames(Mut_forPSI)))]
driver_effect_onPSI=do.call(rbind, driver_effect_onPSI)
driver_effect_onPSI=data.frame(driver_effect_onPSI)
driver_effect_onPSI$fdr=fdrcorr(driver_effect_onPSI$p_value)
driver_effect_onPSI$fdr_corr_pvalue=fdrcorr(driver_effect_onPSI$corr_pvalue)
#Significant genes
paste(unlist(rownames(driver_effect_onPSI)[which((driver_effect_onPSI$fdr_corr_pvalue)<0.1)]), collapse=', ')
#Order by 
driver_effect_onPSI=driver_effect_onPSI[order(abs(1-driver_effect_onPSI$FC), decreasing = T),]
head(driver_effect_onPSI)
# Mutation Profile
na.omit(driver_effect_onPSI[match(APA_mach$GENE.NAME, rownames(driver_effect_onPSI)),])
driver_effect_onPSI[match(APA_mach$GENE.NAME, rownames(driver_effect_onPSI)),]

colnames(driver_effect_onPSI)[2]='FC_PSI(Mutated/WildType)'
cat(rownames(driver_effect_onPSI[which(driver_effect_onPSI$fdr_corr_pvalue<0.1 &
                                   driver_effect_onPSI$`FC_PSI(Mutated/WildType)` >1.2),]), sep=', ')
write.csv(driver_effect_onPSI[,c(2,3,5)],
          '/Users/sinhas8/APA_Adriana/2.Data/drivers_mutation_association_withPSI.csv')
######################
# CORE APA GENES
######################
CORE_APA_mach_mutated=as.numeric(colSums(rowSubset(Mut_forPSI, as.character(CORE_APA_mach$Polyadenylation)), na.rm = T)>0)
err_handle(c(p_value=summary(lm(df2plot_FF$Normalized_PSI ~
                                  CORE_APA_mach_mutated+df2plot_FF$cancer_type))$coefficients[2,4],
             FC=median(df2plot_FF$Normalized_PSI[as.logical(CORE_APA_mach_mutated)], na.rm=T)/
               median(df2plot_FF$Normalized_PSI[!as.logical(CORE_APA_mach_mutated)], na.rm=T))) 

######################
# All APA genes
######################
CORE_APA_mach_mutated=as.numeric(colSums(rowSubset(Mut_forPSI, as.character(APA_mach$GENE.NAME)), na.rm = T)>4)
CORE_APA_mach_mutated=as.numeric(colSums(rowSubset(Mut_forPSI, as.character(CORE_APA_mach$Polyadenylation)), na.rm = T)>2)
err_handle(c(p_value=summary(lm(df2plot_FF$Normalized_PSI ~
                                  CORE_APA_mach_mutated+df2plot_FF$cancer_type))$coefficients[2,4],
             FC=median(df2plot_FF$Normalized_PSI[as.logical(CORE_APA_mach_mutated)], na.rm=T)/
               median(df2plot_FF$Normalized_PSI[!as.logical(CORE_APA_mach_mutated)], na.rm=T))) 
wilcox.test(df2plot$PSI ~ df2plot$race, alternative='l')
# APA_mach_id=na.omit(match(APA_mach$GENE.NAME, rownames(Mut_forPSI)))
# APA_mach_id_onPSI=mclapply(APA_mach_id, function(x) err_handle(c(p_value=wilcox.test(df2plot_FF$Normalized_PSI ~ unlist(Mut_forPSI[x,]))$p.value,
#                                                                        FC=median(df2plot_FF$Normalized_PSI[as.logical(Mut_forPSI[x,])], na.rm=T)/
#                                                                          median(df2plot_FF$Normalized_PSI[!as.logical(Mut_forPSI[x,])], na.rm=T) )),
#                                 mc.cores = detectCores())
# names(APA_mach_id_onPSI)=rownames(Mut_forPSI)[na.omit(match(APA_mach$GENE.NAME, rownames(Mut_forPSI)))]
# APA_mach_id_onPSI=do.call(rbind, APA_mach_id_onPSI)
# APA_mach_id_onPSI=data.frame(APA_mach_id_onPSI)
# APA_mach_id_onPSI$fdr=fdrcorr(APA_mach_id_onPSI$p_value)
# colnames(APA_mach_id_onPSI)[2]='FC_PSI(Mutated/WildType)'
# write.csv(APA_mach_id_onPSI, '/Users/sinhas8/APA_Adriana/2.Data/APA_machgenes_mutation_assocaition_onPSI.csv')

######################
# All APA genes
######################
mach_act_Exp=colMedians(prob$mRNA[APA_mach_id,])
df2plot$APA_machinary_activity=mach_act_Exp[match(df2plot$sample_id, gsub('-','.',prob$samples))]
wilcox.test(df2plot$APA_machinary_activity ~ df2plot$race, alternative='l')
df2plot=na.omit(df2plot)
tiff('/Users/sinhas8/APA_Adriana/APAmachinary_PanCan_Nov16.tif',width = 400, height = 400)
ggplot(df2plot, aes(x=race,y=APA_machinary_activity, fill=race))+
  geom_boxplot(data=df2plot, aes(fill=race))+
  labs(x="Race", y = "APA_machinary_activity")+
  theme_classic(base_size = 20)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)
dev.off()

######################
##Genewise test
######################
Gene_corr_withPSI<-function(geneName){
  mach_act_Exp=prob$mRNA.norm[match(geneName, prob$genes),]
  df2plot$APA_machinary_activity=mach_act_Exp[match(df2plot$sample_id, gsub('-','.',prob$samples))]
  unlist(cor.test(df2plot$PSI, df2plot$APA_machinary_activity, method='spear')[c(3,4)])
}
APA_corr_with_APAgenes=t(sapply(rownames(APA_mach_id_onPSI), function(x) err_handle(Gene_corr_withPSI(x)) ))
write.csv(APA_corr_with_APAgenes,'/Users/sinhas8/APA_Adriana/2.Data/Table0_APA_machgenes_corrwithnPSI.csv')
######################
##Mutation Test
######################
geneList_id=APA_mach_id
mach_act=colSums(prob_wdMut$Mut[geneList_id,])
df2plot$APA_machinary_Mutation=mach_act[match(df2plot$sample_id, gsub('-','.', colnames(prob_wdMut$Mut)))]
fisher.test(cbind(table(as.numeric(as.logical(df2plot$APA_machinary_Mutation[df2plot$race=='AA']))),
                  table(as.numeric(as.logical(df2plot$APA_machinary_Mutation[df2plot$race=='EA'])))  ),
            alternative = 'l')
######################
##FOr our cohort
######################
######################
##APA diff in different manner
######################
PSI_tumor=read.csv('/Users/sinhas8/APA_Adriana/2.Data/PSI_Tumor.csv')
PSI_tumor=PSI_tumor[,-1]
myhead(PSI_tumor)
PSI_Norm=read.csv('/Users/sinhas8/APA_Adriana/2.Data/PSI_Normal.csv')
myhead(PSI_Norm)
PSI_Norm=PSI_Norm[,-1]
delta_PSI=PSI_tumor - PSI_Norm

Demo_Age_sex=read.csv('/Users/sinhas8/APA_Adriana/2.Data/Demo_Age_Sex.csv')
race_df=t(Demo_Age_sex[match(gsub('X','',colnames(delta_PSI)), Demo_Age_sex$PATIENT.ACC.),])[c(8),]
names(race_df)=NULL
demo=read.csv('/Users/sinhas8/APA_Adriana/2.Data/demo.csv')
histology_df=c(t(demo[match(gsub('X','',colnames(delta_PSI)),demo$acc),])[c(2),])
names(histology_df)=NULL
require('robustbase')
delta_PSI_load=colMedians(as.matrix(delta_PSI))
dat=data.frame(delta_PSI=delta_PSI_load,
               Hist=histology_df,
               Race=race_df)
tiff('/Users/sinhas8/APA_Adriana/4.Analysis/NCIMD_PSILoad.tiff')
ggplot(dat, aes(y=delta_PSI, fill=Race, x=Hist))+
  geom_boxplot(outlier.colour = NA)+
  labs(title="Delta-PSI Load NCI-MD",x="Race", y = "Delta-PSILoad")+
  # facet_grid( ~Hist, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4, method.args = list(alternative='l'))+
  guides(fill=FALSE)+
  scale_fill_brewer(palette='Set1')#+
dev.off()

dat=rbind(data.frame(PSI=c(apply(PSI_tumor, 1, function(x) x)),
                     Hist=histology_df,
                     Race=race_df,
                     Sample='Tumor'),
          data.frame(PSI=c(apply(PSI_Norm, 1, function(x) x)),
                     Hist=histology_df,
                     Race=race_df,
                     Sample='Norm'))

dat=na.omit(dat)
dat=dat[dat$Race !='H',]
dat$Race=factor(as.character(dat$Race))
dat=dat[dat$Hist=='adeno' | dat$Hist=='sq',]
dat$Hist=factor(as.character(dat$Hist))
tiff('/Users/sinhas8/APA_Adriana/4.Analysis/NCIMD_PSI_both.tiff')
ggplot(dat, aes(y=PSI, fill=Race, x=Sample))+
  geom_boxplot(outlier.colour = NA)+
  labs(title="Del_PSI NCI-MD",x="Race", y = "Del_PSI")+
  facet_grid( ~ Hist, scales = "free")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette='Set1')
dev.off()
dim(prob$mRNA)