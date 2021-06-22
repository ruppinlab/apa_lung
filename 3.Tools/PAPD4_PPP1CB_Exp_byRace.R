require(ggpubr)
df2plot=rbind(data.frame(Exp=scale(prob$mRNA[which(prob$genes=='PAPD4'),prob$types=='BRCA' & (prob$race=='AA' | prob$race=='EA')]),
                         race=prob$race[prob$types=='BRCA' & (prob$race=='AA' | prob$race=='EA')],
                         geneName='PAPD4'),
              data.frame(Exp=scale(prob$mRNA[which(prob$genes=='PPP1CB'),prob$types=='BRCA' & (prob$race=='AA' | prob$race=='EA')]),
                         race=prob$race[prob$types=='BRCA' & (prob$race=='AA' | prob$race=='EA')],
                         geneName='PPP1CB')
)
df2plot=na.omit(df2plot)

tiff('/Users/sinhas8/APA_Adriana/4.Analysis/PAPD4_PPP1CB_Exp_byRace.tiff')
ggplot(df2plot, aes(x=race, y=Exp, fill=race))+
  geom_boxplot()+
  labs(title="",x="Race", y = "Expression in BRCA")+
  facet_grid( ~geneName, scales = "free", space = "free_x")+
  theme_classic(base_size = 25)+
  stat_compare_means(method = "wilcox.test", label = "p", label.x = 1.5,
                     size = 4)+
  guides(fill=FALSE)+
  scale_fill_brewer(palette='Set1')
dev.off()
