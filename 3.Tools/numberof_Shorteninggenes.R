######################
##APA diff in different manner
######################
source('/Users/sinhas8/myCustom_functions.R')
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
######################
# Previous method
######################
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
######################
# Count the number of number of low/high
######################
dat=rbind(data.frame(PSI=colSums((PSI_tumor)>0.5),
                     Hist=histology_df,
                     Race=race_df,
                     Sample='Tumor'),
          data.frame(PSI=colSums(PSI_Norm>0.5),
                     Hist=histology_df,
                     Race=race_df,
                     Sample='Normal'))

dat=na.omit(dat)
dat=dat[dat$Race !='H',]
dat$Race=factor(as.character(dat$Race))
dat=dat[dat$Hist=='adeno' | dat$Hist=='sq',]
dat$Hist=factor(as.character(dat$Hist))
tiff('/Users/sinhas8/APA_Adriana/4.Analysis/NCIMD_geneswith_3Prime_shortening_bothTUmorNorm_Upd.tiff')
ggplot(dat, aes(y=PSI, fill=Race, x=Sample))+
  geom_boxplot(outlier.colour = NA)+
  labs(title="Delta PSI NCI-MD",x="Race", y = "Genes with 3'-shortening")+
  # facet_grid( ~ Hist, scales = "free")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette='Set1')
dev.off()
aggregate(dat$PSI, list(dat$Race), function(x) median(x))
