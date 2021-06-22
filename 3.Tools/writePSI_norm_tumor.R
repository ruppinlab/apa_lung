setwd('/Users/sinhas8/APA_Adriana/')
mat=read.csv('5.Previous_Versions/Finished_Project_APA_prognostic-Marker/Data/complete_prox_Distal.csv')
mat=mat[,c(1, grep('N',colnames(mat)))]
rownames(mat)=mat[,1]
mat=mat[,-1]
ID_Map=read.csv('2.Data/Id_mapping.csv')
colnames(mat)=ID_Map[!is.na(match(gsub('T','',as.character(unlist(ID_Map[,1]))), gsub('N','',substring(colnames(mat), 2)) )),2]
prox=apply(mat, 1, function(x) sapply(x, function(y) strsplit(y, ';')[[1]][1] ))
dist=apply(mat, 1, function(x) sapply(x, function(y) strsplit(y, ';')[[1]][2] ))
prox=apply(prox, 1, as.numeric)
rownames(prox) = rownames(mat)
dist=apply(dist, 1, as.numeric)
rownames(dist) = rownames(mat)
PSI=(prox + 0.001)/(prox+dist+0.001)
colnames(PSI)=colnames(mat)
write.csv(PSI, '2.Data/PSI_Normal.csv')


# Tumor PSI
mat=read.csv('5.Previous_Versions/Finished_Project_APA_prognostic-Marker/Data/complete_prox_Distal.csv')
mat=mat[,c(1, grep('T',colnames(mat)))]
rownames(mat)=mat[,1]
mat=mat[,-1]
ID_Map=read.csv('2.Data/Id_mapping.csv')
colnames(mat)=ID_Map[!is.na(match(gsub('T','',as.character(unlist(ID_Map[,1]))), gsub('T','',substring(colnames(mat), 2)) )),2]
prox=apply(mat, 1, function(x) sapply(x, function(y) strsplit(y, ';')[[1]][1] ))
dist=apply(mat, 1, function(x) sapply(x, function(y) strsplit(y, ';')[[1]][2] ))
prox=apply(prox, 1, as.numeric)
rownames(prox) = rownames(mat)
dist=apply(dist, 1, as.numeric)
rownames(dist) = rownames(mat)
PSI=(prox + 0.001)/(prox+dist+0.001)
colnames(PSI)=colnames(mat)
write.csv(PSI, '2.Data/PSI_Tumor.csv')

