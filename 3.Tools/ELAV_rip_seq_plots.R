########################################################################
# Introduction
########################################################################
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)
ripReads=import('/Users/sinhas8/Downloads/wgEncodeSunyRipSeqGm12878Elavl1SigRep1.bigWig',format = 'bigwig')
mir_coor=readxl::read_xlsx('/Users/sinhas8/Downloads/mir file for sanju with genomic coordinates (1).xlsx', sheet = 2)
score=ripReads$score
start_coor=ripReads@ranges@start
end_coor=ripReads@ranges@start+ripReads@ranges@start
chr_info=sapply(1:24, function(x) 
  rep(ripReads@seqnames@values[x], ripReads@seqnames@lengths[x]))
ripReads@seqnames@values
chr_info=unlist(chr_info)
chr_info_trimmed=gsub('chr','',chr_info)
########################################################################
# Define func
########################################################################
function_reads<-function(infunc_chr=mir_coor$chr[1],
                         infunc_start=mir_coor$proximal_pos...3[1],
                         infunc_end=mir_coor$distal_pos...5[1]){
  chr_id=chr_info_trimmed==infunc_chr
  range_id=(start_coor<=infunc_start & end_coor<=infunc_start) |
    (start_coor<=infunc_end & end_coor<=infunc_end)
  score[which(chr_id & range_id)]
}
Gm12878_elav=sapply(1:nrow(mir_coor),
       function(x) function_reads(infunc_chr=mir_coor$chr[x],
               infunc_start=mir_coor$proximal_pos[x],
               infunc_end=mir_coor$distal_pos[x]))
ripReads_k562=import('/Users/sinhas8/Downloads/wgEncodeSunyRipSeqK562Elavl1SigRep1.bigWig',
                     format = 'bigwig')
score_k562=ripReads_k562$score
start_coor_k562=ripReads_k562@ranges@start
end_coor_k562=ripReads_k562@ranges@start+ripReads_k562@ranges@start
chr_info_k562=sapply(1:24, function(x) 
  rep(ripReads_k562@seqnames@values[x], ripReads_k562@seqnames@lengths[x]))
chr_info_k562=unlist(chr_info_k562)
chr_info_trimmed_k562=gsub('chr','',chr_info_k562)
function_reads_k562<-function(infunc_chr=mir_coor$chr[1],
                         infunc_start=mir_coor$proximal_pos[1],
                         infunc_end=mir_coor$distal_pos[1]){
  chr_id_k562=chr_info_trimmed_k562==infunc_chr
  range_id_k562=(start_coor_k562<=infunc_start & end_coor_k562<=infunc_start) |
    (start_coor_k562<=infunc_end & end_coor_k562<=infunc_end)
  score_k562[which(chr_id_k562 & range_id_k562)]
}
k562_elav=sapply(1:nrow(mir_coor),function(x)
  function_reads_k562(infunc_chr=mir_coor$chr[x],
                      infunc_start=mir_coor$proximal_pos[x],
                      infunc_end=mir_coor$distal_pos[x]))
########################################################################
# Plot it:: Part 1
########################################################################
combineddf=lapply(1:10, function(x) data.frame(Reads=c(Gm12878_elav[[x]], k562_elav[[x]]),
                                               Sample=c(rep('Tumor', length(Gm12878_elav[[x]])),
                                                        rep('Normal', length(k562_elav[[x]]))),
                                               geneName=mir_coor$gene_name[x]))
combineddf=do.call(rbind, combineddf)
combineddf=combineddf[combineddf$geneName!='MIR99AHG',]
combineddf$geneName=as.character(combineddf$geneName)
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal.tiff', height=600, width=600)
ggplot(combineddf, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  facet_wrap(geneName~., scales = 'free', nrow=3)+
  theme_bw(base_size = 20)+
  stat_compare_means(label = 'p', method='wilcox')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

########################################################################
# Plot it:: Part 2:: # For short in tumor
########################################################################
combineddf_short=lapply(grep('short',mir_coor$gene_class...1), function(x) data.frame(Reads=c(Gm12878_elav[[x]], k562_elav[[x]]),
                                               Sample=c(rep('Tumor', length(Gm12878_elav[[x]])),
                                                        rep('Normal', length(k562_elav[[x]]))),
                                               geneName=mir_coor$gene_name[x]))
combineddf_short=do.call(rbind, combineddf_short)
combineddf_short=combineddf_short[combineddf_short$geneName!='MIR99AHG',]
combineddf_short$geneName=as.character(combineddf_short$geneName)
saveRDS(combineddf_short,'/Users/sinhas8/APA_Adriana/ELAV_mirofInterest.RDS')
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal_Short.tiff', height=300, width=300)
ggplot(combineddf_short, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  theme_bw(base_size = 20)+
  stat_compare_means(label = 'p', method='wilcox')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()
########################################################################
# Plot it:: Part 2:: # For long in tumor
########################################################################
combineddf_long=lapply(grep('long',mir_coor$gene_class...1),
                        function(x) data.frame(Reads=c(Gm12878_elav[[x]], k562_elav[[x]]),
                                               Sample=c(rep('Tumor', length(Gm12878_elav[[x]])),
                                                        rep('Normal', length(k562_elav[[x]]))),
                                               geneName=mir_coor$gene_name[x]))
combineddf_long=do.call(rbind, combineddf_long)
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal_long.tiff',
     height=300, width=300)
ggplot(combineddf_long, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  theme_bw(base_size = 18)+
  stat_compare_means(label = 'p', method='wilcox')+
  annotate(geom="text", x=1.5, y=3,
           label=paste('N/T=',aggregate(combineddf_long$Reads ~ combineddf_long$Sample, combineddf_long, median)[1,2]/
                         aggregate(combineddf_long$Reads ~ combineddf_long$Sample, combineddf_long, median)[2,2]),
           color="red")
dev.off()


########################################################################
# Plot individual miRNA:: miR-155
########################################################################
combineddf_short_MIR155HG=combineddf_short[combineddf_short$geneName %in% 'MIR155HG',]
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal_MIR155HG.tiff',
     height=300, width=300)
ggplot(combineddf_short_MIR155HG, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  facet_wrap(geneName~., scales = 'free', nrow=3)+
  theme_bw(base_size = 20)+
  stat_compare_means(label = 'p', method='wilcox')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()
########################################################################
# Plot individual miRNA:: Let7B
########################################################################
combineddf_short_MIRLET7BHG=combineddf_short[combineddf_short$geneName %in% 'MIRLET7BHG',]
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal_MIRLET7BHG.tiff',
     height=300, width=300)
ggplot(combineddf_short_MIRLET7BHG, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  facet_wrap(geneName~., scales = 'free', nrow=3)+
  theme_bw(base_size = 20)+
  stat_compare_means(label = 'p', method='wilcox')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()
########################################################################
# Plot individual miRNA:: miR-17
########################################################################
combineddf_long_MIR17HG = combineddf[combineddf$geneName %in% 'MIR17HG',]
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal_MIR17HG.tiff',
     height=300, width=300)
ggplot(combineddf_long_MIR17HG, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  facet_wrap(geneName~., scales = 'free', nrow=3)+
  theme_bw(base_size = 20)+
  stat_compare_means(label = 'p', method='wilcox')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

########################################################################
# lincRNAs
########################################################################
linc_coor=readxl::read_xlsx('/Users/sinhas8/Downloads/mir file for sanju with genomic coordinates (1).xlsx', sheet = 1)
linc_coor=linc_coor[which(linc_coor$gene_name %in% c('DLEU1', 'PVT1')),]
colnames(linc_coor)
linc_coor=linc_coor[,c(8, 1, 2, 4)]
x=2
Gm12878_elav=sapply(1:nrow(linc_coor), function(x)
  function_reads(infunc_chr=linc_coor$chr[x],
                 infunc_start=linc_coor$proximal_pos...2[x],
                 infunc_end=linc_coor$distal_pos...4[x]))
k562_elav=sapply(1:nrow(linc_coor), function(x)
  function_reads_k562(infunc_chr=linc_coor$chr[x],
                      infunc_start=linc_coor$proximal_pos...2[x],
                      infunc_end=linc_coor$distal_pos...4[x]))
combineddf_linc=lapply(1:2, function(x) data.frame(Reads=c(Gm12878_elav[[x]], k562_elav[[x]]),
                                               Sample=c(rep('Tumor', length(Gm12878_elav[[x]])),
                                                        rep('Normal', length(k562_elav[[x]]))),
                                               geneName=linc_coor$gene_name[x]))
combineddf_linc=do.call(rbind, combineddf_linc)
head(combineddf_linc)
table(combineddf_linc$geneName)
# 
combineddf_linc_PVT1 = combineddf_linc[combineddf_linc$geneName %in% 'PVT1',]
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal_PVT1.tiff',
     height=300, width=300)
ggplot(combineddf_linc_PVT1, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  facet_wrap(geneName~., scales = 'free', nrow=3)+
  theme_bw(base_size = 20)+
  stat_compare_means(label = 'p', method='wilcox')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

# 
combineddf_linc_DLEU1 = combineddf_linc[combineddf_linc$geneName %in% 'DLEU1',]
tiff('/Users/sinhas8/APA_Adriana/ELAVReads_tumorvsNormal_DLEU1.tiff',
     height=300, width=300)
ggplot(combineddf_linc_DLEU1, aes(y=log(Reads, 10), x=Sample, fill=Sample))+
  geom_boxplot()+
  facet_wrap(geneName~., scales = 'free', nrow=3)+
  theme_bw(base_size = 20)+
  stat_compare_means(label = 'p', method='wilcox')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()
