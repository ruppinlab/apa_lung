##########################################################################################
# TCGA_proteomics_for_BRCA
##########################################################################################
require(statar)
brca_prot=get(load('/Users/sinhas8/Downloads/tcga_brca_prot.RData'))
[grep('ACE2',brca_prot$geneid$symbol
load('/Users/sinhas8/ISLE-Basic/data/TCGA.RData')
mRNA_matched=prob$mRNA[,match(brca_prot$pheno$patient_id, tolower(prob$samples))]
mRNA_matched_SF=mRNA_matched[match(brca_prot$geneid$symbol, prob$genes),]
x=1
df2plot=read.csv('/Users/sinhas8/Downloads/2.Data/TCGA_PSI.csv')
df2plot_matched=df2plot[match(brca_prot$pheno$patient_id, gsub('\\.','-',tolower(df2plot$sample_id))),]
K=5
df2plot_matched$PSI_qunatized=xtile(df2plot_matched$PSI, K)
# When PSI is low - Longer 3'UTR
corr_matrix=sapply(1:nrow(mRNA_matched_SF), function(x) err_handle(unlist(cor.test(mRNA_matched_SF[x,df2plot_matched$PSI_qunatized==1],
                                                                                   brca_prot$expr[x,df2plot_matched$PSI_qunatized==1])[c(3, 4)])) )
corr_matrix=do.call(rbind, corr_matrix)
corr_matrix=data.frame(corr_matrix)
median(corr_matrix$estimate.cor, na.rm = T)

# PSI is high - Shorter 3'UTR
corr_matrix=sapply(1:nrow(mRNA_matched_SF), function(x) err_handle(unlist(cor.test(mRNA_matched_SF[x,df2plot_matched$PSI_qunatized==K],
                                                                                   brca_prot$expr[x,df2plot_matched$PSI_qunatized==K])[c(3, 4)])) )
corr_matrix=do.call(rbind, corr_matrix)
corr_matrix=data.frame(corr_matrix)
median(corr_matrix$estimate.cor, na.rm = T)
# Compute the above in a geneWise manner
K=1000
geneWise_PSI=data.frame(gene=apa_dist[[3]][[6]],
                        PSI_med=rowMeans(apa_dist[[3]][[1]], na.rm = T) )
geneWise_PSI_ranked=geneWise_PSI[order(geneWise_PSI$PSI_med),]
Lowest_PSI_genes= head(geneWise_PSI_ranked$gene, K)
Highest_PSI_genes= tail(geneWise_PSI_ranked$gene, K)
corr_matrix_lowPSI=sapply(match(Lowest_PSI_genes, brca_prot$geneid$symbol),
                   function(x) err_handle(unlist(cor.test(mRNA_matched_SF[x,],
                                                          brca_prot$expr[x,])[c(3, 4)])) )
corr_matrix_lowPSI=do.call(rbind, corr_matrix_lowPSI)
rownames(corr_matrix_lowPSI)=Lowest_PSI_genes
corr_strength_lowPSI=unlist(data.frame(corr_matrix_lowPSI)[,2])

corr_matrix_highPSI=sapply(match(Highest_PSI_genes, brca_prot$geneid$symbol),
                   function(x) err_handle(unlist(cor.test(mRNA_matched_SF[x,],
                                                          brca_prot$expr[x,])[c(3, 4)])) )
corr_matrix_highPSI=do.call(rbind, corr_matrix_highPSI)
corr_strength_highPSI=unlist(data.frame(corr_matrix_highPSI)[,2])

df2plot1=rbind(data.frame(corr_strength=corr_strength_highPSI, PSI_type='High PSI', GeneName=Highest_PSI_genes),
              data.frame(corr_strength=corr_strength_lowPSI, PSI_type='Low PSI', GeneName=Lowest_PSI_genes))
df2plot1$cohort= 'TCGA BRCA'

##########################################################################################
# CCLE proteomics
##########################################################################################
setwd('/Users/sinhas8/APA_Adriana/')
ccle_corr=readxl::read_xlsx('CCLEcorr_protein_mRNA.xlsx', sheet = 2)
mean_PDUI_ccle=readxl::read_xlsx('CCLE_meanPDUI_Supp_Tables.xlsx', sheet = 2)
mean_PDUI_ccle$mean_CCLE_PSI=1-mean_PDUI_ccle$Tumor_mean_PDUI

# Sort by mean PSI
mean_PDUI_ccle=mean_PDUI_ccle[order(mean_PDUI_ccle$mean_CCLE_PSI, decreasing = T),]
head(mean_PDUI_ccle)
df2plot=rbind(data.frame(mRNAvsProt_Cor=ccle_corr$Spearman[
  na.omit(match(unique(mean_PDUI_ccle$Gene_Symbol[mean_PDUI_ccle$Trends=='Shortening']), ccle_corr$`Gene Symbol`))],
  APA='Shortening',
  GeneName=ccle_corr$`Gene Symbol`[
    na.omit(match(unique(mean_PDUI_ccle$Gene_Symbol[mean_PDUI_ccle$Trends=='Shortening']), ccle_corr$`Gene Symbol`))],
  cohort='CCLE'),
  data.frame(mRNAvsProt_Cor=ccle_corr$Spearman[
    na.omit(match(unique(mean_PDUI_ccle$Gene_Symbol[mean_PDUI_ccle$Trends=='Lengthening']), ccle_corr$`Gene Symbol`))],
  APA='Lengthening',
  GeneName=ccle_corr$`Gene Symbol`[
    na.omit(match(unique(mean_PDUI_ccle$Gene_Symbol[mean_PDUI_ccle$Trends=='Lengthening']), ccle_corr$`Gene Symbol`))],
  cohort='CCLE')
  )
colnames(df2plot1)=colnames(df2plot)
levels(df2plot$APA)=c('High PSI', 'Low PSI')
df2plot_comb=rbind(df2plot, df2plot1)
write.csv(df2plot_comb,'/Users/sinhas8/APA_Adriana/2.Data/SuppTable12and13apaContext_mRNAvsProtein_corr.csv')

tiff('/Users/sinhas8/APA_Adriana/ProtvsmRNA_givenAPAcontext.tiff', width = 600, height = 400)
ggplot(df2plot_comb, aes(x=APA, y=mRNAvsProt_Cor, fill=APA))+
  geom_boxplot()+
  facet_wrap(~cohort)+
  stat_compare_means(method.args = list(alternative = "less"))+
  theme_bw(base_size = 18)+
  labs(x='APA Type', y='Correlation Strength\nmRNA vs Protein')+
  theme(axis.text.x =  element_blank())
dev.off()  
write.csv(df2plot, '/Users/sinhas8/APA_Adriana/ProtvsmRNA_givenAPAcontext.csv')
