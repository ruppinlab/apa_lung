# Plot ELAV1 binding
setwd('/Users/sinhas8/Downloads/')
reads=read.csv('wgEncodeSunyRipSeqGm12878Elavl1Pk.broadPeak', sep='\t', header = F)
Mir_of_interest=readxl::read_xlsx('mir file for sanju with genomic coordinates.xlsx', sheet = 2)
tail(Mir_of_interest)
reads$V1=substring(reads$V1, 4,)
Mir_of_interest=Mir_of_interest[,c(2,3,5)]
Mir_of_interest$chr=paste('chr',Mir_of_interest$chr, sep='')
head(Mir_of_interest)

library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

try1=import('/Users/sinhas8/Downloads/wgEncodeSunyRipSeqGm12878Elavl1SigRep1.bigWig',
            format = 'BigWig')

sapply(1:1000, function(x) try1[x,1]@seqnames@values)
mir_reads_mapping=lapply(1:nrow(Mir_of_interest),function(x) 
  intersect(match(Mir_of_interest$chr[x], try1@seqinfo@seqnames),
       c(which(try1@ranges@start< 
               Mir_of_interest$proximal_pos[x] &
               try1@ranges@start+try1@ranges@width <
               Mir_of_interest$proximal_pos[x]),
       which(try1@ranges@start< 
               Mir_of_interest$distal_pos[x] &
               try1@ranges@start+try1@ranges@width <
               Mir_of_interest$distal_pos[x]))))

mir_reads_count=lapply(mir_reads_mapping, function(x) try1$score[x])
sapply(mir_reads_count, sum)