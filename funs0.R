numlassoCalFun <- function(zX,Y, family0,cF0,flag,ll1){

  for (cc0 in ll1:1){
    kk0 <- numlasso2(zX,Y,family0,cF0*flag+cc0-1,colM,nId)
    numS1 <- kk0$num
    if (flag>=numS1[1]){
      rrNum <- kk0$num
      idGM <- kk0$id
      modelNTM<-kk0$mNTM
      modelMis<-kk0$mMis
      betaE <- kk0$beta_est
      break
    }
  }
  if (cc0==1){
    rrNum <- kk0$num
    idGM <- kk0$id
    modelNTM<-kk0$mNTM
    modelMis<-kk0$mMis
    betaE <- kk0$beta_est
  }
  return(list(num=rrNum,id=idGM, mNTM=modelNTM,mMis=modelMis,beta_est = betaE))
}

regfun_fast<- function(xj,Y,family0){
  Z_R <- ifelse(xj==2,1,0);
  Z_A <- xj
  Z_D <- ifelse(xj==0,0,1);
  gR <-summary(fastglm(matrix(cbind(1,Z_R),ncol=2),Y,family=family0))$coefficients
  pR <- gR[2,4]
  gA <-summary(fastglm(matrix(cbind(1,Z_A),ncol=2),Y,family=family0))$coefficients
  pA <- gA[2,4]
  gD <-summary(fastglm(matrix(cbind(1,Z_D),ncol=2),Y,family=family0))$coefficients
  pD <- gD[2,4]
  p_min <- min(c(pR,pA,pD))
  out = list(dosageP = pA, p_min3 = p_min)
  out
} 


GMscreenFun2 <-function(dataX,Y,colM,flag){
  ddX <- dataX[,colM]
  p2 <-dim(ddX)[2]
  zX <- matrix(rep(0,p2*2*dim(ddX)[1]),ncol=2*p2)
  for (j in 1:p2){
    #  zX[,(j-1)*3+1]=0.5*dat[,j]; # additive
    zX[,(j-1)*2+1]=ifelse(ddX[,j]==2,1,0);# recessive
    zX[,(j-1)*2+2]=ifelse(ddX[,j]==0,0,1);# dominant
  }
  
  # ss1 <- numlasso(zX,Y,family0,flag,colM)[1]
  ss20 <- numlasso2(zX,Y,family0,2*flag,colM,nId)
  if (flag==0){
    out1 = ss20
  }else{
    ss2 <- ss20$num[1]
    ss30 <- numlasso2(zX,Y,family0,floor(1.5*flag),colM,nId)
    ss3 <- ss30$num[1]

    if (flag>ss2){
      ll1 <- length((2*flag):(3*flag));
      out1 <- numlassoCalFun(zX,Y, family0,cF0=2,flag,ll1)
    } else {
      if (flag<ss3){
        ll1 <- length((flag+1):floor(1.5*flag));
        out1 <- numlassoCalFun(zX,Y, family0,cF0=1,flag,ll1)
      }else{
        ll1 <- length((floor(1.5*flag)+1):(2*flag));
        out1 <- numlassoCalFun(zX,Y, family0,cF0=1.5,flag,ll1)
      }
    }
  }
  return(out1)
}



createData <- function(dd0,n,cM,nId,nperc){
  # nperc: the percentage of training data
  #cM: the number of features included 
  
  rowId <- sample(1:dim(dd0)[1],n)
  columnId <- sample(1:dim(dd0)[2],cM);
  resX <- sampleCausal(dd0[rowId,columnId],0.01,nId);
  X01 <- resX$cauData
  
  columnId <- columnId[-resX$cauId]
  data00 <- dd0[rowId,columnId]
  dataX <- cbind(X01,data00);
  
  
  if(nId==10){
    Z_R1 <- ifelse(dataX[,5]==2,1,0);
    Z_R2 <- ifelse(dataX[,6]==2,1,0);
    Z_R3 <- ifelse(dataX[,7]==2,1,0);
    Z_D1 <- ifelse(dataX[,8]==0,0,1);
    Z_D2 <- ifelse(dataX[,9]==0,0,1);
    Z_D3 <- ifelse(dataX[,10]==0,0,1);
    XX <- cbind(dataX[,1:4],Z_R1,Z_R2,Z_R3,Z_D1,Z_D2,Z_D3);
  }else if(nId==20){
    Z_R1 <- ifelse(dataX[,11]==2,1,0);
    Z_R2 <- ifelse(dataX[,12]==2,1,0);
    Z_R3 <- ifelse(dataX[,13]==2,1,0);
    Z_R4 <- ifelse(dataX[,14]==2,1,0);
    Z_R5 <- ifelse(dataX[,15]==2,1,0);
    Z_D1 <- ifelse(dataX[,16]==0,0,1);
    Z_D2 <- ifelse(dataX[,17]==0,0,1);
    Z_D3 <- ifelse(dataX[,18]==0,0,1);
    Z_D4 <- ifelse(dataX[,19]==0,0,1);
    Z_D5 <- ifelse(dataX[,20]==0,0,1);
    XX <- cbind(dataX[,1:10],Z_R1,Z_R2,Z_R3,Z_R4,Z_R5,Z_D1,Z_D2,Z_D3,Z_D4,Z_D5);
  }else{
    Z_R1 <- ifelse(dataX[,1]==2,1,0);
    Z_R2 <- ifelse(dataX[,2]==2,1,0);
    Z_R3 <- ifelse(dataX[,3]==2,1,0);
    Z_D1 <- ifelse(dataX[,4]==0,0,1);
    Z_D2 <- ifelse(dataX[,5]==0,0,1);
    Z_D3 <- ifelse(dataX[,6]==0,0,1);
    XX <- cbind(Z_R1,Z_R2,Z_R3,Z_D1,Z_D2,Z_D3,dataX[,7:nId]);
  }
  
  
  Y <- as.matrix(XX)%*%matrix(c(rep(1,nId)),ncol=1)+rnorm(dim(XX)[1],0,1)
  if (nperc ==1){
    return(list(X=dataX, Y=Y)) 
  }else{
    n1 <- floor(n*nperc)
    return(list(X.train=dataX[1:n1,], Y.train=Y[1:n1], X.test = dataX[(n1+1):n,],Y.test=Y[(n1+1):n])) 
  }
}



transLasso<-function(zX,Y,family0,num,colM){
  res11 <- lassoID(zX,Y,1:dim(zX)[1],1:dim(zX)[2],family0,num)
  
  rId <- numeric();
  for (kk in 1:length(colM)){
    rId[2*(kk-1)+1]=2*(colM[kk]-1)+1
    rId[2*(kk-1)+2]=2*(colM[kk]-1)+2
  }
  
  pp <- ifelse(rId[res11[1,]]%%2==0,rId[res11[1,]]%/%2,rId[res11[1,]]%/%2+1)
  pp1 <- unique(pp)
  pp1
}

numlasso2<- function(zX,Y,family0,num,colM,nId){
  #10(4add 3rec 3dom), 20 (10add, 5 rec, 5dom), others (3rec,3dom,else add) 
  if(nId==10){
    modN <- c(4,7) 
  }else if (nId==20){
    modN <- c(10,15) 
  }else{
    modN <-'error'
    break
  }
    
  
  res11 <- lassoID(zX,Y,1:dim(zX)[1],1:dim(zX)[2],family0,num)
  
  rId <- numeric();
  for (kk in 1:length(colM)){
    rId[2*(kk-1)+1]=2*(colM[kk]-1)+1
    rId[2*(kk-1)+2]=2*(colM[kk]-1)+2
  }
  
  pp <- ifelse(rId[res11[1,]]%%2==0,rId[res11[1,]]%/%2,rId[res11[1,]]%/%2+1)
  pp1 <- unique(pp)
  
  Idd <-rId[res11[1,]]
  beta_est <- rep(0,2*nId)     
  
  
  modelNTM <-rep(0,3)
  modelMis <-rep(0,9)
  for (kk in 1:nId){
    yy1 = which(Idd==2*(kk-1)+1)
    yy2 = which(Idd==2*(kk-1)+2)
 
    beta_est[2*(kk-1)+1] = ifelse(length(yy1)==0, 0, res11[2,yy1])
    beta_est[2*(kk-1)+2] = ifelse(length(yy2)==0, 0, res11[2,yy2])
 
    
    if (kk<=modN[1]){
      if (length(yy1)>0&length(yy2)>0){
        modelNTM[1]=modelNTM[1]+1
      } else if(length(yy1)>0&length(yy2)==0){
        modelMis[1]=modelMis[1]+1;#REC
      } else if(length(yy1)==0&length(yy2)>0){
        modelMis[2]=modelMis[2]+1; #DOM
      } else{
        modelMis[3]=modelMis[3]+1; #Miss
      }
    }
    if (kk>modN[1]&kk<=modN[2]){ #REC
      if (length(yy1)>0&length(yy2)==0){
        modelNTM[2]=modelNTM[2]+1
      } else if (length(yy1)>0&length(yy2)>0){
        modelMis[4]=modelMis[4]+1; #ADD
      } else if (length(yy1)==0&length(yy2)>0){
        modelMis[5]=modelMis[5]+1; #DOM
      } else {
        modelMis[6]=modelMis[6]+1; #Miss
      }
    }
    if (kk>modN[2]){ 
      if (length(yy1)==0&length(yy2)>0){
        modelNTM[3]=modelNTM[3]+1 #DOM
      } else if (length(yy1)>0&length(yy2)>0){
        modelMis[7]=modelMis[7]+1 #ADD
      } else if (length(yy1)>0&length(yy2)==0){
        modelMis[8]=modelMis[8]+1 #REC
      } else {
        modelMis[9]=modelMis[9]+1 #Miss
      }
    }
  } 
  
  
  a01 <- intersect(causalId,pp1);
  
  numS0 <- c(length(pp1),length(a01));
  
  res1 <- list(num=numS0,id=a01, mNTM=modelNTM, mMis=modelMis, beta_est = beta_est)
  res1
}




modeDetect2 <- function(selectId,case){
  # case 1: 1-4(ADD),5-7(REC),8-10(DOM)
  # case 2: 
  if (case==1){
    numREC <- sum(selectId==5)+sum(selectId==6)+sum(selectId==7)
    numDOM <- sum(selectId==8)+sum(selectId==9)+sum(selectId==10)
    numADD <- sum(selectId==1)+sum(selectId==2)+sum(selectId==3)+sum(selectId==4) 
    rate <- c(numADD,numREC,numDOM) 
  }else{
    numREC <- sum(selectId==11)+sum(selectId==12)+sum(selectId==13)+sum(selectId==14)+sum(selectId==15)
    numDOM <- sum(selectId==16)+sum(selectId==17)+sum(selectId==18)+sum(selectId==19)+sum(selectId==20)
    numADD =0
    for (i1 in 1:10){
      numADD <- numADD+sum(selectId==i1)
    }
    rate <- c(numADD,numREC,numDOM) 
  }
  rate
}



sampleCausal <- function(dat,r,num){
  # r: minor allele frequency
  # num: the number of causal SNPs selected
  nCol <- dim(dat)[2]
  nrow <- dim(dat)[1]
  columnId <- 1:nCol
  new <-numeric();
  cId <- numeric();
  
  for (i in 1:num){
  ratioR <- 0;
  while (ratioR <r){
    j0 <- sample(columnId,1)
    r0 <- sum(dat[,j0]==0)/nrow
    r1 <- sum(dat[,j0]==1)/nrow
    r2 <- sum(dat[,j0]==2)/nrow
    ratioR <- min(r0,r1,r2)
    if (ratioR<r){
      columnId <- columnId[-j0]
    } else {
      aa=dat[,j0]
    
      if (r0<r2){
        aa=ifelse(aa==2,3,aa);
        aa=ifelse(aa==0,2,aa);
        aa=ifelse(aa==3,0,aa);
        r0 <- sum(aa==0)/nrow
        if (r0<r1){
          aa=ifelse(aa==1,3,aa);
          aa=ifelse(aa==0,1,aa);
          aa=ifelse(aa==3,0,aa);
        }
      }
    }
  }
  new <- cbind(new,aa)
  cId <- c(cId,j0)
  }
  L <- list(cauData=new, cauId=cId)
  L
}

lassoID <- function(data,Y,rowId,columnId,family0,flag){
  # data: SNPs on a chromnsome
  # family : binomial , gaussian
  # flag: 0(lambda.1se); num(the number of features desired to be selected)
  
  
  X <-  as.matrix(data[rowId,columnId])
  y <- Y[rowId]
  gla <- cv.glmnet(X, y, family=family0,nfolds = 10)
  # fit <- glmnet(X,y,family=family0,alpha=alpha0,lambda=gla$lambda.1se)
  #  Id <-columnId[which(fit$beta!=0)]
  # coeff <- fit$beta[which(fit$beta!=0)]
  if (flag==0){
    cc <- gla$lambda.1se
  } else {
    if (flag>=max(gla$nzero)){
      cc <- gla$lambda[length(gla$lambda)]
    }else{
      lambdaId = which(flag<gla$nzero)[1]-1
      cc <- gla$lambda[lambdaId]
    }
  }

  coef.beta = coef(gla, s = cc)
  # coef.beta = coef(gla, s = "lambda.1se") 
  Id <-columnId[which(coef.beta[2:length(coef.beta)]!=0)]
  a <- coef.beta[2:length(coef.beta)]
  coeff <- a[which(coef.beta[2:length(coef.beta)]!=0)]
  
  res <- rbind(Id,coeff)
  
  return(res)  
}

