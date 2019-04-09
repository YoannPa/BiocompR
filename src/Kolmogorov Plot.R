##Kolmogorov Plot

##Imports
library(doParallel)
library(broom)
library(reshape2)
library(corrplot)

#Set WD
setwd("/media/yoann/3CAD87DD271F7BEC/PhD_2018-2021/Hematopoiesis_DNMT1_Hypomorph/data/")

#Set Number of cores to be used
registerDoParallel(cores=3)

##Load matrix of beta values from CpGs or DMRs
Beta.Mat<-t(readRDS("HSPC_DMRs_NoSNP_Cov5.RDS"))

#Give the type of regions ("DMRs" or "CpGs")
region_type<-"DMRs"

#Are the samples Merged? (TRUE or FALSE)
is_merged<-F

#Create list of sample beta values
list_samples <- lapply(seq_len(nrow(Beta.Mat)), function(i) Beta.Mat[i,])
names(list_samples)<-rownames(Beta.Mat)

#Remove mpp34 samples: #After discussions: Do not keep mpp34 samples
if(is_merged == F) {
  list_samples$mpp34_1<-NULL
  list_samples$mpp34_2<-NULL
  list_samples$mpp34_3<-NULL
} else {
  list_samples$mpp34<-NULL
}

##Make all pairwise combinations
table_combinations<-expand.grid(names(list_samples),names(list_samples))

##Kolmogorov-Smirnov 2 samples test of distribution (if all CpGs and 3 threads time ~ 2 hours)
List_ks.tests<-foreach(comb=seq(nrow(table_combinations))) %dopar%
  ks.test(list_samples[[table_combinations[comb,1]]],
          list_samples[[table_combinations[comb,2]]])

#Create Table
List_ks.tests<-lapply(List_ks.tests,tidy)
Table.ks.tests<-as.data.frame(do.call("rbind",List_ks.tests))[,c(1,2)]

#Create Table of P-Value
Table.ks.pval<-cbind(table_combinations,Table.ks.tests$p.value)
#Create Table of Maximum Deviations
Table.ks.D<-cbind(table_combinations,Table.ks.tests$statistic)

#Cast a molten tables into matrices
Mat.ks.pval<-dcast(Table.ks.pval,Var1 ~ Var2)
rownames(Mat.ks.pval)<-Mat.ks.pval$Var1
Mat.ks.pval<-as.matrix(Mat.ks.pval[,-1])

Mat.ks.D<-dcast(Table.ks.D,Var1 ~ Var2)
rownames(Mat.ks.D)<-Mat.ks.D$Var1
Mat.ks.D<-as.matrix(Mat.ks.D[,-1])

if(is_merged == F) {
  Mat.ks.D<-Mat.ks.D[c(1:5,8,6,7,9:11,14,12,13,15:17,20,18,19,21:23,26,24,25,27:36),
                     c(1:5,8,6,7,9:11,14,12,13,15:17,20,18,19,21:23,26,24,25,27:36)]
  
  Mat.ks.pval<-Mat.ks.pval[c(1:5,8,6,7,9:11,14,12,13,15:17,20,18,19,21:23,26,24,25,27:36),
                           c(1:5,8,6,7,9:11,14,12,13,15:17,20,18,19,21:23,26,24,25,27:36)]
  } else {
  Mat.ks.D<-Mat.ks.D[c(1,3,2,4,6,5,7,9,8,10,12,11,13:16),
                     c(1,3,2,4,6,5,7,9,8,10,12,11,13:16)]
  
  Mat.ks.pval<-Mat.ks.pval[c(1,3,2,4,6,5,7,9,8,10,12,11,13:16),
                     c(1,3,2,4,6,5,7,9,8,10,12,11,13:16)]
}

#Adjust colors with quantiles
quantile(Table.ks.tests$statistic, seq(0, 1, 0.1))

#Plot the Kolmogorov Deviations ################################################
corrplot(Mat.ks.D, title = "Pairwise Kolmogorov Cumulative Relative Frequencies Maximum Deviations",
         method = "color", # addCoef.col = "white",
         tl.col = "black", tl.srt = 45,
         diag = FALSE,
         order = "original",
         col=c(rep("white",100),
               rep("#045a8d",3),rep("#0570b0",1),rep("#3690c0",2),
               rep("#74a9cf",1),rep("#fdd49e",1),
               rep("#fdbb84",2),rep("#fc8d59",9),rep("#ef6548",3),
               rep("#d7301f",3),rep("#7f0000",6),rep("#7c0303",71),rep("#6a3d9a",1)),
         cl.pos = "r", cl.ratio = 0.1, cl.lim = c(0,0.31),
         tl.offset = 0.4, addgrid.col = "grey",
         mar=c(0,0,2,0))

## Correlation plot with Associated Standard Errors ############################
corrplot(round(Mat.ks.pval,2),
         title = "Pairwise Kolmogorov-Smirnov Test P-Values",
         method = "color",type = "lower",addCoef.col = "black",tl.col = "black",
         col = c(rep("black",100),rep("lavenderblush",1),rep("#d1e5f0",4),
                 rep("#a6d96a",5),rep("#60bb5d",82),rep("#66bd63",8)),
         diag = FALSE, order = "original", number.cex= 0.6,
         cl.lim = c(min(Mat.ks.pval),1), cl.ratio = 0.1,cl.pos = "b", mar=c(0,0,2,0))







#Ggplot the Kolmogorov Deviations ##############################################






