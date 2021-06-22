# Load LUAD protein-mRNA correlation
genes_APA=read_xlsx('/Users/sinhas8/Downloads/Supp Tables Jan 27 (1).xlsx',
                    sheet=4, col_names = T)
corr_tumor=read_xlsx('/Users/sinhas8/APA_Adriana/protein_mRNA_correlation_tumorvsNormal.xlsx',
                     sheet = 1)

sum(genes_APA$...25=='repressed')*0.1

shortened_genes=head(genes_APA$...5[genes_APA$...25=='repressed'], 
                     25)
elongated_genes=tail(genes_APA$...5[genes_APA$...25=='enhanced'],
                     25)
df2plot=rbind(data.frame(cor=corr_tumor$cor[corr_tumor$Gene %in% shortened_genes],
                      Type='High PSI'),
           data.frame(cor=corr_tumor$cor[corr_tumor$Gene %in% elongated_genes],
                      Type='Low PSI'))
df2plot$cohort='LUAD TCGA'
tiff('/Users/sinhas8/APA_Adriana/ProtvsmRNA_givenAPAcontext_Lung_Tumor.tiff',
     height = 300, width = 300)
ggplot(df2plot, aes(x=Type, fill=Type, y=cor))+
  geom_boxplot()+
  facet_wrap(cohort~.)+
  stat_compare_means(method = 'wilcox')+
  labs(y='Correlation Strength between\n mRNA and Protein levels')+
  theme_bw(base_size = 15)
dev.off()

