setwd("path0") # please put your path here
source('funs.R')
#######################################################################################
#Simulate GWAS dataset using GWAsimulator.exe (provided by C. Li @vanderbilt U)       #
# usage: in the file GWAsimulator_v2.0_windows, press Shift and right click the mouse,# 
# then choose 'open command window here'                                              #
#Winthin the command window, type GWAsimulator control.dat 1, where 1 is a random seed#
#######################################################################################

#------------------------Simulation procedure-----------------------------------------#
library(glmnet)
library(LDheatmap)
library(combinat)
library(parallel)

cl <- makeCluster(5)
clusterEvalQ(cl,library(fastglm))


# generate simulation dataset, 
d1 <-read.table('chr1.dat');



set.seed(0)
dd0 <- clean(d1,0.0009)
# clean data, drop those SNP with minor allele frequency <0.03; 
# made the coding as follows: 2: AA, 1:Aa, 0:aa. A is the minor allele frequency



M=1000;    
n <- 3000;
#n <-5000;
cM <-dim(dd0)[2]; # the number of features
family0 = "gaussian"
p0 = 0.01
nId = 20;
causalId <- 1:nId
#causalId <- 1:6;
d=ceiling(n/log(n));

train_perc = 0.8


set.seed(0)
res <- matrix(rep(0,3*M),ncol=3)
FDR <- matrix(rep(0,3*M),ncol=3)
numMode <- matrix(rep(0,9*M),ncol=9)
modeMis <- matrix(rep(0,9*M),ncol=9)
modeNTM <- matrix(rep(0,3*M),ncol=3)
betaE_GM <- matrix(rep(0,M*2*nId),nrow=M)
betaE_LA <- matrix(rep(0,M*nId),nrow=M)
betaE_CATT <- matrix(rep(0,M*nId),nrow=M)

d1 <- ceiling((1-train_perc)*n/log((1-train_perc)*n))
FPR <- 0.0005



for (i in 1:M){
  data_i= createData(dd0,n,cM,nId,nperc=train_perc)
  dataX = data_i$X.train
  Y = data_i$Y.train
  
  # train step

  res1 <- parSapply(cl,dataX, FUN=regfun_fast, Y=Y, family0=family0)
  pM <- as.numeric(res1['p_min3',])
  
  
  pA <- as.numeric(res1['dosageP',])
  pv <- sort(pA,decreasing = FALSE, index.return = TRUE)
  L1 <- which(pv$x>p0)[1]
  
  flag <- min(c(cM*FPR+nId,d1,L1))
  
  
  psort<- sort(pM, decreasing = FALSE, index.return = TRUE)
  Idall <- psort$ix[psort$x<p0]
  if (length(Idall)>d){
    Idall <- Idall[1:d]
  }
  colM <- Idall # select the top d markers
  

  dataX.test = data_i$X.test
  Y.test = data_i$Y.test
  
  
  # GMscreen
  res_GM <- GMscreenFun2(dataX.test,Y.test,colM,ceiling(flag))
  rrNum <- res_GM$num
  idGM <- res_GM$id
  modelNTM<- res_GM$mNTM
  modelMis<- res_GM$mMis
  betaE_GM[i,] <- res_GM$beta_est
  
  
  
  
  # Lasso
  
  Idall <- pv$ix[1:d]
  ddata <- dataX.test[,Idall];
  res22 <- lassoID(ddata,Y.test,1:dim(ddata)[1],1:dim(ddata)[2],family0,flag)
  a22 <- intersect(causalId,Idall[res22[1,]])
  
  beta_la <- rep(0,nId)
  
  Idall2 <- Idall[res22[1,]]
  for (kk in 1:nId){
    yy1 = which(Idall2==kk)
    
    beta_la[kk] = ifelse(length(yy1)==0, 0, res22[2,yy1])
  }
  betaE_LA[i,] = beta_la 
  
  #  CATT
  d1 = max(c(dim(res22)[2],rrNum[1]))
  d2 = ifelse(flag==0,d1,flag)
  aV <- pv$ix[1:d2]
  a11 <- intersect(causalId,aV)  
  
  datX2 = cbind(Y.test, dataX.test[,aV])
  colnames(datX2) <-c('Y1', as.character(aV))
  f1 = summary(glm(Y1~.,data=datX2,family=family0))
  CC1 = f1$coefficients[2:dim(f1$coefficients)[1],]
  
  beta_CA = rep(0,nId)
  for(j2 in 1:nId){
    id0 = which(aV==j2)
    beta_CA[j2] =ifelse(length(id0)==0,0,ifelse(id0>dim(CC1)[1],0,CC1[id0,'Estimate']))
  }
  betaE_CATT[i,] = beta_CA 
  
  modeNum3 <- modeDetect2(a22,2)
  modeNum2 <- modeDetect2(a11,2)
  modeNum1 <- modeDetect2(idGM,2)
  
  res[i,]=c(rrNum[2],length(a22),length(a11))
  FDR[i,]=c(rrNum[1]-rrNum[2],dim(res22)[2]-length(a22),length(aV)-length(a11))
  numMode[i,]=c(modeNum1,modeNum3,modeNum2)
  modeMis[i,]=modelMis
  modeNTM[i,]=modelNTM
}
