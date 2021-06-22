require(data.table)
pdui_GTEX=fread('/Users/sinhas8/Downloads/PDUI.txt')
pdui_GTEX=data.frame(pdui_GTEX)

geneInfo=pdui_GTEX[,1:3]
demo_gtex=get(load('/Users/sinhas8/Downloads/pheno.RData'))
demo_gtex=demo_gtex[match(sapply(colnames(pdui_GTEX), function(x) paste0(strsplit(x, '\\.')[[1]][1:2], collapse = '-') ),
                          demo_gtex$SUBJID),]
pdui_GTEX_FF=pdui_GTEX[,which(demo_gtex$RACE==2 | demo_gtex$RACE==3)]
demo_gtex=demo_gtex[which(demo_gtex$RACE==2 | demo_gtex$RACE==3),]
pdui_GTEX_FF_scaled=apply(pdui_GTEX_FF, 1, scale)
pdui_GTEX_FF_scaled=t(pdui_GTEX_FF_scaled)
PDUI_load=colMedians(pdui_GTEX_FF_scaled, na.rm = T)

wilcox.test(PDUI_load~ factor(demo_gtex$RACE), alternative='l')
