##APA effect survival
##Libraries needed
require(survival)
#Loading files
setwd('/Users/sinhas8/APA_Adriana/')
mat=read.csv('2.Data/AZ_T_N_2.complete_tab.csv')
colnames(mat)
##Matrix Preprocessing
rownames(mat)=mat[,1]
mat=mat[,-1]
ID_Map=read.csv('2.Data/Id_mapping.csv')
colnames(mat)=ID_Map[!is.na(match(as.character(unlist(ID_Map[,1])), substring(colnames(mat), 2))),2]
prox=apply(mat, 1, function(x) sapply(x, function(y) strsplit(y, ';')[[1]][1] ))
dist=apply(mat, 1, function(x) sapply(x, function(y) strsplit(y, ';')[[1]][2] ))

prox=apply(prox, 1, as.numeric)
rownames(prox) = rownames(mat)
dist=apply(dist, 1, as.numeric)
rownames(dist) = rownames(mat)


PSI=(prox + 0.001)/(prox+dist+0.001)

colnames(PSI)=colnames(mat)
#png('/home/sinhas/APA_dist.png')
#plot(density(PSI, na.rm=T))
#dev.off()
#write.csv(PSI,'/home/sinhas/PSI.csv')
####Cox regression 
dd=read.csv('/Users/sinhas8/APA_Adriana/2.Data/demo.csv')
dd=dd[na.omit(match(colnames(PSI), dd$acc)),]
PSI_filtered_for_cox=PSI#[,-ncol(PSI)]

######Demographics:: 
surv=dd[,c(7, 10)]
hist=dd$hist
stage=dd$stage

Demo_Age_sex=read.csv('/Users/sinhas8/APA_Adriana/2.Data/Demo_Age_Sex.csv')
gender = Demo_Age_sex$gender
age    = Demo_Age_sex$age
race   = Demo_Age_sex$race
rownames(PSI_filtered_for_cox)= sapply(rownames(PSI_filtered_for_cox), function(x) strsplit(x, ' ')[[1]][1])


APA_on_Survival<-function(mat=PSI_filtered_for_cox, Candidate_Genes=rownames(PSI_filtered_for_cox), surv=surv, stage=stage, hist=hist){
	##APA_on_Survival
	cox_results=apply(mat[Candidate_Genes,], 1, function(x) coxph(Surv(surv$survival, surv$lungcancer_death_all_years) ~ x + hist + stage + gender + age + race))
	cox_summary=sapply(cox_results, function(x) summary(x)$coefficients['x',c(1,5)])
	cox_summary=t(cox_summary)
	cox_summary=cox_summary[order(cox_summary[,2]),]
	FDR=p.adjust(cox_summary[,2], method='fdr')
	cbind(cox_summary, FDR)
}

#temp=APA_on_Survival(mat=PSI_filtered_for_cox, Candidate_Genes=rownames(PSI_filtered_for_cox), surv=surv, stage=stage, hist=hist)
#write.csv(temp, '/home/sinhas/APA_on_Survival_after_hist_stage_age_sex_race.csv')

#Top_1k_Variable= order(apply(PSI_filtered_for_cox, 1, function(x) var(x, na.rm=T)), decreasing=T)[1:1000]
#temp1=APA_on_Survival(mat=PSI_filtered_for_cox, Candidate_Genes=Top_1k_Variable, surv=surv, stage=stage, hist=hist)
#write.csv(temp1, '/home/sinhas/APA_on_Survival_Top_1k_Variable_Genes_after_hist_stage_age_sex_race.csv')


#Variance_order= order(apply(PSI_filtered_for_cox, 1, function(x) var(x, na.rm=T)), decreasing=T)
#write.csv(PSI[Variance_order,],'/home/sinhas/PSI_ordered_by_Variance.csv')


##########Correlational Matrix
#png('/home/sinhas/All_genes_APA_corr-matrix.png')
#corrplot(PSI_Corr, diag = FALSE, order = "hclust", tl.pos = "td", tl.cex = 0.2, method = "color", type = "upper")
#dev.off()


##Preprocessing
#Genelist=read.csv('/home/sinhas/APA_Adriana/2.Data/new_genes.csv')
#Genelist=as.character(unlist(Genelist))
#PSI_filtered_for_cox=PSI_filtered_for_cox[,]

#Risk Groups
M=50
K=98

require(survriskpred)
projectPath <- '/home/sinhas/temp_delete/'
outputName <- "SurvivalRiskPrediction"

exprTrain=PSI_filtered_for_cox[,-(M:K)]
exprTest=as.data.frame(PSI_filtered_for_cox[,M:K])

covTrain=data.frame(age=age[-(M:K)])
covTest=data.frame(age=age[M:K])
geneId=data.frame(Symbol=as.character(unlist((rownames(exprTrain)))))

tme=surv[,2]
status=surv[,1]
#Function
resList2 <- survRiskPredict(exprTrain, covTrain, exprTest,covTest,
                           geneId, status, tme, geneSelect ="pc",
                           nriskgroups = 2, progIndexPerc = 50,
                           cvMethod = "10fold", nperm = 0, landmarktime = 0,
                           alpha = .001, ncomp = 2, mixing = 1, pcrgenes = 10,
                           projectPath = projectPath,
                           outputName = outputName)

############################################################################################################################
####SuperPC
############################################################################################################################
data=list(x=as.matrix(exprTrain), y=tme[-(M:K)], censoring.status= status[-(M:K)], featurenames=geneId)
train.obj<- superpc.train(data, type="survival")
cv.obj<-superpc.cv(train.obj, data)

superpc.plotcv(cv.obj)

#Thresh=0.9
Thresh=0.7
data.test=list(x=as.matrix(exprTest), y=tme[(M:K)], censoring.status= status[(M:K)], featurenames=geneId)
lrtest.obj<-superpc.lrtest.curv(train.obj, data,data.test)
superpc.plot.lrtest(lrtest.obj)
fit.cts<- superpc.predict(train.obj, data, data.test, threshold=Thresh, n.components=1, prediction.type="continuous")
superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.pred)

png('/home/sinhas/KP_Train_New_Genes.png')
##Kaplan Mier of Training data
fit.groups<- superpc.predict(train.obj, data, data, threshold=Thresh, n.components=2, prediction.type="discrete")
superpc.fit.to.outcome(train.obj, data, fit.groups$v.pred)
fit<-survfit(Surv(data$y,data$censoring.status)~fit.groups$v.pred)
plot(fit, col=2:3, xlab="time", ylab="Prob survival")
dev.off()

##Kaplan Mier of Testing set
#png('/home/sinhas/KP_Test.png')

fit.groups<- superpc.predict(train.obj, data, data.test, threshold=Thresh, n.components=1, prediction.type="discrete")
superpc.fit.to.outcome(train.obj, data.test, fit.groups$v.pred)
plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.pred), col=2:3, xlab="time", ylab="Prob survival")
#dev.off()

png('/home/sinhas/test1B.png')
fit.red<- superpc.predict.red(train.obj, data, data.test, threshold=Thresh)
fit.redcv<- superpc.predict.red.cv(fit.red, cv.obj,  data,  threshold=Thresh)
superpc.plotred.lrtest(fit.redcv)
write.csv(superpc.listfeatures(data, train.obj, fit.red, num.features=10), '/home/sinhas/test1.csv')
dev.off()




##Kaplan Mier
##TRAINING
fit.groups<- superpc.predict(train.obj, data, data, threshold=Thresh, n.components=1, prediction.type="discrete")
SurvTrain=data.frame(time=data$y, status=data$censoring.status, group=fit.groups$v.pred)
superpc.fit.to.outcome(train.obj, data, fit.groups$v.pred)
#fit<-survfit(Surv(data$y,data$censoring.status)~fit.groups$v.pred)
fit<-survfit(Surv(time,status)~discrete, data=SurvTrain)
ggsurvplot(
  fit, 
  data=SurvTrain,
  size = 1,                 # change line size
  palette = 
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
#  risk.table = TRUE,        # Add risk table
# risk.table.col = "strata",# Risk table color by groups
  legend.labs = 
    c("Subtype I", "Subtype II"),    # Change legend labels
#  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
ggsave('/home/sinhas/KP_Train_New_genes.png')

##TESTING
fit.groups<- superpc.predict(train.obj, data, data.test, threshold=Thresh, n.components=1, prediction.type="discrete")
SurvTest=data.frame(time=data.test$y, status=data.test$censoring.status, group=fit.groups$v.pred)
#fit<-survfit(Surv(data$y,data$censoring.status)~fit.groups$v.pred)
fit<-survfit(Surv(time,status)~discrete, data=SurvTest)
ggsurvplot(
  fit, 
  data=SurvTest,
  size = 1,                 # change line size
# palette =     c("#E7B800", "#2E9FDF"),# custom color palettes
# conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
# risk.table = TRUE,        # Add risk table
# risk.table.col = "strata",# Risk table color by groups
  legend.labs =  c("Subtype I", "Subtype II"),    # Change legend labels
# risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

ggsave('/home/sinhas/KP_Test_New_Genes.png')
