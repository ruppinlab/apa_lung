##Our cohort-level PSI
source('/Users/sinhas8/myCustom_functions.R')
PSI_tumor=read.csv('/Users/sinhas8/APA_Adriana/5.Previous_Versions/Finished_Project_APA_prognostic-Marker/PSI_new.csv')
row.names(PSI_tumor)=(PSI_tumor)[,1]
PSI_tumor=(PSI_tumor)[,-1]
Demo_Age_sex=read.csv('/Users/sinhas8/APA_Adriana/2.Data/Demo_Age_Sex.csv')
race_df=t(Demo_Age_sex[match(gsub('X','',colnames(PSI_tumor)), Demo_Age_sex$PATIENT.ACC.),])[c(8),]
names(race_df)=NULL
demo=read.csv('/Users/sinhas8/APA_Adriana/2.Data/demo.csv')
histology_df=c(t(demo[match(gsub('X','',colnames(PSI_tumor)),demo$acc),])[c(2),])
names(histology_df)=NULL

df2plot=data.frame(PSI=colMeans(PSI_tumor, na.rm = T),
                             race=race_df,
                             cancer_type=histology_df)

df2plot=df2plot[df2plot$race !='H',]
df2plot$Race=factor(as.character(df2plot$race))
df2plot=df2plot[df2plot$cancer_type=='adeno' | df2plot$cancer_type=='sq',]
df2plot$cancer_type=factor(as.character(df2plot$cancer_type))
levels(df2plot$cancer_type)=c('LUAD', 'LUSC')
write.csv(df2plot, '2.Data/OurCohort_PSI.csv')
######################
##Plto Figures
######################
tiff('/Users/sinhas8/APA_Adriana/PSI_byCancerType_OurCohort.tif',
     width = 600, height = 400)
ggplot(df2plot, aes(x=as.character(cancer_type),y=PSI, fill=race))+
  geom_boxplot(data=df2plot, aes(fill=race))+
  labs(title="PSI in our Cohort",x="Race", y = "PSI")+
  facet_grid( ~cancer_type, scales = "free", space = "free_x")+
  theme_bw(base_size = 20)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4, method.args = list(alternative='l'))+
  guides(fill=FALSE)+
  scale_fill_brewer(palette=colorType_Set)
dev.off()
