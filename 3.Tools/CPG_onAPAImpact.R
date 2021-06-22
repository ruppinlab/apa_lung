######Step 1: Find CpGs per gene
require(data.table)
require(statar)
require(ggplot2)
setwd('/Users/sinhas8/Downloads/')
gene_pos=read.csv('mart_export.txt', sep = '\t')
cpg=read.csv('CpG_genome', sep='\t', header=F)
head(cpg)
count_cpgs<-function(start=gene_pos$Gene.start..bp., end=gene_pos$Gene.end..bp.){
  #sum(between(cpg$V2, start-min(1000, start), end) | between(cpg$V3, start-min(1000, start), end) )
  sum(between(cpg$V2, start, end) | between(cpg$V3, start, end) )
}
dim(gene_pos)
numCpG=apply(gene_pos, 1, function(x) count_cpgs(as.numeric(x[4]), as.numeric(x[5])) )
length=apply(gene_pos, 1, function(x) as.numeric(x[5]) -as.numeric(x[4]) )
gene_wdCpGs=data.frame(geneName=gene_pos$Gene.name, numCpG, length)
dim(gene_wdCpGs)
#####Step 2: Find CpGs per gene
corr_wdAPA=read.csv('Correlation_and_pvalues.csv')
gene_wdCpGs=gene_wdCpGs[na.omit(match(corr_wdAPA$gene_name, gene_wdCpGs$geneName)),]
corr_wdAPA=corr_wdAPA[!is.na(match(corr_wdAPA$gene_name, gene_wdCpGs$geneName)),]
df2plot=data.frame(corr=corr_wdAPA$spr_corr, CpGs=gene_wdCpGs$numCpG, 
                   CpGs_distinct=factor(xtile(gene_wdCpGs$numCpG, 20)),
                   len=gene_wdCpGs$length, 
                   len_distinct=factor(xtile(gene_wdCpGs$length, 10)))
p1=ggplot(df2plot, aes(y=corr, x=len_distinct))+geom_boxplot()+
  labs(x='Length of Gene', y='Impact (Corr RHo) on Expression')
p2=ggplot(df2plot, aes(y=corr, x=CpGs_distinct))+geom_boxplot()+
  labs(x='Number of CpGs', y='Impact (Corr RHo) on Expression')

tiff('/Users/sinhas8/APA_Adriana/APAImpactExp_corr_wdgeneLenght.tiff')
grid.arrange(p1, p2, ncol=2)
dev.off()
###
cor.test(corr_wdAPA$spr_corr, gene_wdCpGs$numCpG)
TukeyHSD(aov(df2plot$corr ~  df2plot$CpGs_distinct))
cor.test(df2plot$corr, df2plot$CpGs, method='s' )
TukeyHSD(aov(df2plot$corr ~  df2plot$len_distinct))
cor.test(df2plot$corr, gene_wdCpGs$length, method='s' )
summary(lm(df2plot$corr ~  df2plot$CpGs_distinct+df2plot$len_distinct))
