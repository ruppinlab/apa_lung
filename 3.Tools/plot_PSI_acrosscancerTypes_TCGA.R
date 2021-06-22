######################
##APA across race in TCGA cancer types
######################
##APA Breast cancer
source('/Users/sinhas8/myCustom_functions.R')
require(parallel); require(rowr);require(ggplot2); require(ggpubr)
colorType_Set='Set1'
genelist=read.csv('/Users/sinhas8/Downloads/Correlation_and_pvalues.csv')
setwd('/Users/sinhas8/Downloads/')
file_list=list.files()[grep('Combined_PDUIs.txt$',list.files())]
length(file_list)
prob=readRDS('/Users/sinhas8/Downloads/TCGA_withMut.RDS')
prob$race=factor(prob$race)
levels(prob$race)=c('AI', 'AS', 'AA', 'AS', 'HN', 'EA')

more_shortening_in_EA <- function(filename_infunc=file_list[23]){
  cancer_apa=read.csv(filename_infunc, sep='\t')
  # cancer_apa[,1]=sapply(cancer_apa[,1], function(x) strsplit(as.character(x), '\\|')[[1]][2])
  # cancer_apa=cancer_apa[na.omit(match(genelist$gene_name, cancer_apa[,1])),]
  geneNames=sapply(strsplit(as.character(cancer_apa[,1]), '\\|'), function(x) x[2])
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
  return(list(cancer_apa, infunc_race, infunc_sex, infunc_age, infunc_stage, geneNames,
              rep(filename_infunc, length(infunc_stage))) )
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

# systematic charting of APA events in TCGA
# including names
apa_dist=mclapply(1:length(file_list), 
                  function(x) err_handle(more_shortening_in_EA(file_list[x])), mc.cores = detectCores())
apa_mat_TCGA=do.call(cbind.fill, lapply(apa_dist, function(x) x[[1]] ))
race_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[2]]))) ))
sex_TCGA=unlist(lapply(apa_dist, function(x) as.character(unlist(err_handle(x[[3]]))) ))
age_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[4]])) ))
stage_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[5]])) ))
geneNames=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[6]])) ))
cancertype_TCGA=unlist(lapply(apa_dist, function(x) unlist(err_handle(x[[7]])) ))

df2plot=data.frame(sample_id=colnames(apa_mat_TCGA),
                   PSI=colMeans(1-apa_mat_TCGA, na.rm = T),
                   race=race_TCGA,
                   sex=sex_TCGA,
                   age=age_TCGA,
                   stage=stage_TCGA,
                   cancer_type=cancertype_TCGA)
# df2plot=na.omit(df2plot)
# df2plot$cancer_type=sapply(df2plot$cancer_type, function(x) strsplit(as.character(x), '_')[[1]][2] )
tcga=read.csv('/Users/sinhas8/Project_Chromotrypsis/Results_New/Nov_28/Supp_Table3.csv')
df2plot$info_tcga.Tissue_Type=tcga$info_tcga.Tissue_Type[match(df2plot$cancer_type, tcga$info_tcga.hist)]
df2plot$info_tcga.CellofOrigin=tcga$info_tcga.CellofOrigin[match(df2plot$cancer_type, tcga$info_tcga.hist)]
write.csv(df2plot, '2.Data/TCGA_PSI.csv')
######################
##Plto Figures
######################
df2plot=read.csv('/Users/sinhas8/Downloads/2.Data/TCGA_PSI.csv')
df2plot=na.omit(df2plot)
colnames(df2plot)
levels(df2plot$cancer_type)
tiff('/Users/sinhas8/APA_Adriana/PSI_byCancerType_Nov15.tif',width = 1800, height = 650)
ggplot(df2plot, aes(x=cancer_type,y=PSI, fill=race))+
  geom_boxplot(data=df2plot, aes(fill=race))+
  labs(title="PSI in 23 cancer types (TCGA)",x="Race", y = "PSI")+
  facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)
dev.off()

# Figure 4B: Pan Cancer Figure
tiff('/Users/sinhas8/APA_Adriana/PSI_byRace_PanTCGA.tif',width = 400, height = 400)
ggplot(df2plot, aes(x=race,y=PSI, fill=race))+
  geom_boxplot(data=df2plot, aes(fill=race))+
  labs(title="PSI in PanCan (TCGA)",x="Race", y = "PSI")+
  # facet_grid( ~info_tcga.Tissue_Type + info_tcga.CellofOrigin, scales = "free", space = "free_x")+
  theme_classic(base_size = 20)+
  # stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
  #                    size = 4)+
  annotate(geom = 'text', x = 1.5, y = 0.5,label=
             paste('P=', summary(lm(df2plot$PSI ~ df2plot$race+df2plot$sex+df2plot$age+df2plot$info_tcga.Tissue_Type))$coefficients[2,4], sep='')) +
  # guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)
dev.off()
