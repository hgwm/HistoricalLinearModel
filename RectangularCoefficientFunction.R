rm(list = ls())

#import functions
library(glmnet)
library(Matrix)
library(ggplot2)
library(lubridate)
library(caTools)
library(forecast)
library(zoo)
library(ranger)
library(missRanger)
library(GPfit)
library(lattice)
library(RColorBrewer)
library(segmented)

get_index_mat <- function(M = 4,max_lag=3) {
  K <- (M+1)*max_lag # total number of nodes
  what<-expand.grid(0:(max_lag-1),0:M)
  data.frame("Node" = 1:K, "Y_coord" = what$Var1, "X_coord" = what$Var2)
}


get_D <- function(M,max_lag) {
  ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  K <- nrow(ind_mat)
  num_pairs <- (M+1)*(max_lag-1)
  
  # First we build the D_h matrix
  D_v <- Matrix(0, nrow = num_pairs, ncol = K)
  ind_diff_v <- which(diff(ind_mat[,3]) == 0)
  
  D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v])] <- -1
  D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v]+1)] <- 1
  D_v<-rbind(D_v,sparseMatrix(i=1:(M+1),j=which(ind_mat$Y_coord == max(ind_mat$Y_coord)), x=rep(1,(M+1)),dims = c(M+1,K) )) # add vanishing boundary to top
  
  # Secondly we have the D_v matrix
  num_pairs <- (M)*(max_lag)
  D_h <- Matrix(0, nrow = num_pairs, ncol = K)
  ind_mat2 <- ind_mat[order(ind_mat[,2]),]
  ind_diff_h<- which(diff(ind_mat2[,2]) == 0)
  D_h[cbind(1:num_pairs, sort(ind_mat2[ind_diff_h,1]))] <- -1
  D_h[cbind(1:num_pairs,sort(ind_mat2[ind_diff_h+1,1]))] <- 1
  D_h<-rbind(D_h,sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == min(ind_mat$X_coord)), x=rep(1,max_lag),dims = c(max_lag,K) ) + sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == max(ind_mat$X_coord)), x=rep(-1,max_lag),dims = c(max_lag,K) )) # add periodicity condition
  
  return(list(D_h,D_v))
}


get_D2 <- function(M,max_lag) {
  ind_mat <- get_index_mat(M = M,max_lag = max_lag)
  K <- nrow(ind_mat)
  
  
  # First we build the D_h matrix
  ind_diff_v2 <- which(diff(diff(ind_mat[,3])) == 0)
  num_pairs <- length(ind_diff_v2)
  
  D_v <- Matrix(0, nrow = num_pairs, ncol = K)
  if(length(ind_diff_v2)>0){
    D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v2])] <- 1
    D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v2]+1)] <- -2
    D_v[cbind(1:num_pairs, ind_mat$Node[ind_diff_v2]+2)] <- 1
  }
  D_v<-rbind(D_v,sparseMatrix(i=1:(M+1),j=which(ind_mat$Y_coord == max(ind_mat$Y_coord)), x=rep(1,(M+1)),dims = c(M+1,K) )) # add vanishing boundary to top
  
  # Secondly we have the D_v matrix
  ind_mat2 <- ind_mat[order(ind_mat[,2]),]
  ind_diff_h2<- which(diff(diff(ind_mat2[,2])) == 0)
  num_pairs <- length(ind_diff_h2)
  D_h <- Matrix(0, nrow = num_pairs, ncol = K)
  D_h[cbind(1:num_pairs, sort(ind_mat2[ind_diff_h2,1]))] <- 1
  D_h[cbind(1:num_pairs,sort(ind_mat2[ind_diff_h2+1,1]))] <- -2
  D_h[cbind(1:num_pairs,sort(ind_mat2[ind_diff_h2+2,1]))] <- 1
  D_h<-rbind(D_h,sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == min(ind_mat$X_coord)), x=rep(1,max_lag),dims = c(max_lag,K) ) + sparseMatrix(i=1:max_lag,j=which(ind_mat$X_coord == max(ind_mat$X_coord)), x=rep(-1,max_lag),dims = c(max_lag,K) )) # add periodicity condition
  
  return(list(D_h,D_v))
}


get_grp_normsH <- function(bvec, M,max_lag) {
  grps <- as.factor(rep(0:M, rep(max_lag,M+1)))
  norms <- tapply(bvec, grps, FUN = function(x){
    if(length(x)>1){
      #cummax(abs(as.numeric(runquantile(x[length(x):1],k=5,endrule = "quantile",probs = 0.5,align = "left"))))[length(x):1]
      #cumsum(abs(as.numeric(x[length(x):1])))[length(x):1]
      #cumsum(as.numeric(runmean(x,k=5)[length(x):1])^2)[length(x):1]
      cumsum(as.numeric(x[length(x):1])^2)[length(x):1]
    } else{
      x
    }
  })
  
  as.numeric(do.call("c", norms))
}


plot_bvec_better<-function(bvec,M,max_lag){
  
  #plot ground truth beta
  index_dat<-get_index_mat(M,max_lag)
  delta<-compute_delta(bvec,M,max_lag)
  
  mat<-matrix(0,ncol = M+1,nrow = max_lag)
  mat[cbind(index_dat$Y_coord+1,index_dat$X_coord+1)]<-bvec
  mat<-mat[1:min(ceiling(max(delta)/10)*10,M+1),]
  mat<-mat[,c((M+2-max_lag):(M+1),1:(M-max_lag+1))]
  
  coul <- colorRampPalette(brewer.pal(8, "Blues"))(30)
  coul2 <- colorRampPalette(brewer.pal(8, "Reds"))(8)
  colorvec<-c(coul2[1],"#FFFFFF",coul[3:length(coul)])
  levelplot(t(mat),col.regions = colorvec,at=sort(c(min(min(mat),-min(abs(mat[mat!=0])*2)),-1e-10,seq(from=1e-10,to=max(max(mat),1e-5),length.out=length(colorvec)-1 ))),
            ylab=list("Lag (Days)",cex=1.4),aspect = "xy",xlab=list("Day of the year" ,cex=1.4),colorkey=list(axis.text=list(cex=1.2)),
            scales=list(x=list(cex=1.2),y=list(cex=1.2)))
}


get_x<-function(x_ts,M,max_lag,first_day){
  K=(M+1)*max_lag
  nrowX<-(length(x_ts)-max_lag+1)
  i_vec<-rep(0, max_lag*nrowX)
  j_vec<-rep(0, max_lag*nrowX)
  x_vec<-rep(0, max_lag*nrowX)
  
  indexCounter<-0
  for(i in 1:nrowX){
    i_vec[(indexCounter+1):(indexCounter+max_lag)]<-rep(i,max_lag)
    j_vec[(indexCounter+1):(indexCounter+max_lag)]<-(((indexCounter + (first_day-1)*max_lag) %% K):((indexCounter+max_lag-1 + (first_day-1)*max_lag) %% K))+1
    x_vec[(indexCounter+1):(indexCounter+max_lag)]<-x_ts[(i+max_lag-1):(i)]
    indexCounter<-indexCounter+max_lag
  }
  
  sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(nrowX, K))
}

compute_delta<-function(bvec,M,max_lag){
  deltavec<-rep(0,M+1)
  
  index_dat<-get_index_mat(M,max_lag)
  
  for(i in 0:M){
    curobs<-bvec[index_dat$X_coord==i]
    curobs<-curobs[length(curobs):1]
    deltavec[i+1]<-max_lag-max(which(cumsum(curobs)==0),0)
  }
  deltavec
}



optimize_hp<-function(hp_range,XTX,hTh,vTv,pTp,XTy,Xval,ytrain,yval,n_hp,To_Zero){
  n_init<-30
  R2<-rep(0,round(n_init))
  N_train<-sum(!To_Zero)/length(ytrain)
  #N_train<-length(ytrain)/sum(!To_Zero)
  #N_train<-1
  
  max_w_h<-hp_range[[1]][2]
  max_w_v<-hp_range[[2]][2]
  
  min_w_h<-hp_range[[1]][1]
  min_w_v<-hp_range[[2]][1]
  
  w_h<-runif(length(R2),min=min_w_h,max=max_w_h)
  w_v<-runif(length(R2),min=min_w_v,max=max_w_v)
  
  
  for(h in 1:length(R2)){
    pred_beta<-rep(0,K)
    pred_beta[!To_Zero] <- Matrix::solve(XTX+ exp(w_h[h])*N_train*hTh + exp(w_v[h])*N_train*vTv, XTy, sparse = TRUE)
    Y_hat_val=Xval %*% pred_beta
    R2[h]<-1-sum((yval- Y_hat_val)^2)/sum((mean(yval)-yval)^2)
    print(paste0("This is R2[h]: ",R2[h], " , and h is: ",h," and w_h is: ",w_h[h], " and w_v is: ", w_v[h]))
  }
  
  while(h<n_hp){
    
    w_h_scaled<-(w_h- min_w_h)/(max_w_h- min_w_h)
    w_v_scaled<-(w_v- min_w_v)/(max_w_v- min_w_v)

    mod<-GP_fit(X=data.frame(w_h=w_h_scaled, w_v=w_v_scaled),Y=R2)
    
    w_h_long<-runif(1000,min=min_w_h,max=max_w_h)
    w_v_long<-runif(1000,min=min_w_v,max=max_w_v)
    
    w_h_long_scaled<-(w_h_long- min_w_h)/(max_w_h- min_w_h)
    w_v_long_scaled<-(w_v_long- min_w_v)/(max_w_v- min_w_v)
    
    pred <- predict.GP(mod, xnew = data.frame(w_h=w_h_long_scaled,w_v=w_v_long_scaled))
    mu <- pred$Y_hat
    sigma <- sqrt(pred$MSE)
    Z <- (mu - max(pred$Y_hat))/sigma
    expected_imp <- sigma*(Z  * pnorm(Z) + dnorm(Z))
    expected_imp[is.na(expected_imp)]<-0
    
    nextone<-which.max(expected_imp)
    w_h<-c(w_h,w_h_long[nextone])
    w_v<-c(w_v,w_v_long[nextone])
    
    h=h+1
    pred_beta<-rep(0,K)
    pred_beta[!To_Zero] <- Matrix::solve(XTX+ exp(w_h[h])*N_train*hTh + exp(w_v[h])*N_train*vTv, XTy, sparse = TRUE)
    Y_hat_val=Xval %*% pred_beta
    R2<-c(R2,1-sum((yval- Y_hat_val)^2)/sum((mean(yval)-yval)^2))
    print(paste0("Gaussian process: This is R2[h]: ",R2[h], " , and h is: ",h," and w_h is: ",w_h[h], " and w_v is: ", w_v[h]))
  }
  list(w_h,w_v,R2)
}

get_elbow<-function(qVec,R2Vec){
  qVecScaled<-(qVec-min(qVec))/(max(qVec)-min(qVec))
  R2VecScaled<- (R2Vec-min(R2Vec))/(max(R2Vec)-min(R2Vec))
  qseq<-seq(from=0,to=1,length.out=1000000)
  R2Line<-1- qseq
  distF0Scaled<-rep(0,length(R2Vec))
  for(i in 1:length(distF0Scaled)){
    distF0Scaled[i]<- min((qVecScaled[i]-qseq)^2 + (R2VecScaled[i]-R2Line)^2) * sign(qVecScaled[i]-qseq)[which.min((qVecScaled[i]-qseq)^2 + (R2VecScaled[i]-R2Line)^2)]
  }
  distF0Scaled
}

get_elbow_menger<-function(qVec,R2Vec){
  qVecscaled<-(qVec-min(qVec))/(max(qVec)-min(qVec))
  R2Vecscaled<- (R2Vec-min(R2Vec))/(max(R2Vec)-min(R2Vec))
  R2Vecscaled<- R2Vecscaled[order(qVecscaled)]
  qVecscaled<-sort(qVecscaled)
  distF0Scaled<-rep(0,length(R2Vec))
  firstq<-which.min(qVecscaled)
  lastq<-which.max(qVecscaled)
  for(i in 1:length(distF0Scaled)){
    triangMat<-matrix(c(qVecscaled[c(firstq,i,lastq)],R2Vecscaled[c(firstq,i,lastq)],rep(1,3)),nrow=3)
    bot<-(dist(triangMat[c(1,2),])*dist(triangMat[c(1,3),])*dist(triangMat[c(2,3),]))
    if(bot==0){
      distF0Scaled[i]<-0
    } else{
      distF0Scaled[i]<- 2*abs(det(triangMat))/bot
    }
  }
  #distF0Scaled<-runmean(distF0Scaled,k=17) #smoothen output
  qVecSorted<-sort(qVec)
  R2VecSorted<-R2Vec[order(qVec)]
  list(qVecSorted[qVecSorted>=0.1 & qVecSorted<= 0.9],R2VecSorted[qVecSorted>=0.1 & qVecSorted<= 0.9],distF0Scaled[qVecSorted>=0.1 & qVecSorted<= 0.9])
}




compute_beta_R2<-function(est,true){
  1-sum((est-true)^2)/sum((true-mean(true))^2)
}

optimize_q<-function(grpnorms,Nq=500,Xmat_train,Xmat_val,y_train,y_val,hMat,vMat,w_h_best,w_v_best){
  q_ngn<-seq(from=0,to=1-365/dim(Xmat_train)[2],length.out=100) #important to select this carfully
  R2<-rep(0,length(q_ngn))
  for(h in 1:length(R2)){
    To_Zero<-grpnorms<quantile(grpnorms,q_ngn[h])
    p_over_n<-sum(!To_Zero)/length(y_train)
    #p_over_n<-length(y_train)/sum(!To_Zero)
    #p_over_n<-1
    XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
    hProd_sparse<-crossprod(hMat[,!To_Zero])
    vProd_sparse<-crossprod(vMat[,!To_Zero])
    
    pred_beta<-rep(0,length(pred_beta))
    pred_beta[!To_Zero]<- Matrix::solve(XXProd_sparse+ w_h_best*p_over_n*hProd_sparse + w_v_best*p_over_n*vProd_sparse, XyProd_sparse, sparse = TRUE)
    Y_hat_val=Xmat_val %*% pred_beta
    R2[h]<-1-sum((y_val- Y_hat_val)^2)/sum((y_val-mean(y_val))^2)
    print(paste0("This is sparse R2[h]: ",R2[h], " , and h is: ",h))
  }
  
  distance2<-get_elbow_menger(q_ngn,R2)
  whichbestq_menger<-which.max(distance2[[3]])
  finalq<-distance2[[1]][whichbestq_menger]
  
  searchRange<-c(max(0,finalq-0.05),min(1,finalq+0.05))
  
  q_ngn2<-runif(Nq- length(R2),min=searchRange[1],max=searchRange[2]) #important to select this carfully
  R22<-rep(0,length(q_ngn2))
  for(h in 1:length(R22)){
    To_Zero<-grpnorms<quantile(grpnorms,q_ngn2[h])
    p_over_n<-sum(!To_Zero)/length(y_train)
    #p_over_n<-length(y_train)/sum(!To_Zero)
    #p_over_n<-1
    XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
    XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
    hProd_sparse<-crossprod(hMat[,!To_Zero])
    vProd_sparse<-crossprod(vMat[,!To_Zero])
    
    pred_beta<-rep(0,length(pred_beta))
    pred_beta[!To_Zero]<- Matrix::solve(XXProd_sparse+ w_h_best*p_over_n*hProd_sparse + w_v_best*p_over_n*vProd_sparse, XyProd_sparse, sparse = TRUE)
    Y_hat_val=Xmat_val %*% pred_beta
    R22[h]<-1-sum((y_val- Y_hat_val)^2)/sum((y_val-mean(y_val))^2)
    if(h%%50==0) print(h/length(R22))
  }
  q_all<-c(q_ngn,q_ngn2)
  R2_all<-c(R2,R22)
  
  return(list(q_all,R2_all))
}


###########################################################################################################
######################################   preprocess data   ###################################################
############################################################################################################


gridcodes<-c(1344, 1405, 1409, 134, 303, 320, 355, 2052, 2042, 1795, 1759)

#get rainfall data
DailyTotalprcp <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/ClimateData/ClimateData/DailyTotalprcp.csv")
Daily_tmean <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/ClimateData/ClimateData/Daily_tmean.csv")

#gridcodes<-as.numeric(DailyTotalprcp[,1])


for(i in 1:length(gridcodes)){
  gridcode<-gridcodes[i]
  
  Precip<-as.numeric(DailyTotalprcp[DailyTotalprcp$Allresults==gridcode,2:ncol(DailyTotalprcp)])
  Temp<-as.numeric(Daily_tmean[Daily_tmean$Allresults==gridcode,2:ncol(Daily_tmean)])
  Rain<-Precip
  Rain[Temp<=0]<-0
  Rain[Rain<0.1]<-0
  
  #toadd<-min(Rain[Rain!=0],na.rm = T)/2
  #Rain[Rain==0]<-toadd
  #Rain<-log(Rain)
  
  #remove leap days from Rain
  #Date_seq<-seq(from=as.Date("1979-10-01"),to=as.Date("2018-09-30"),by="day")
  #Rain<-Rain[seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day") %in% Date_seq]
  Date_seq<-seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day")
  leapdays<-which(day(Date_seq)==29 & month(Date_seq)==2)
  before_lds<-leapdays-1
  Rain[before_lds]<-(Rain[before_lds]+Rain[leapdays])/2
  Rain<-Rain[-leapdays]
  #Rain[Rain<0.1]<-0 #ask ali about this
  
  #put into matrix form
  simu_x<-matrix(Rain,byrow = T,ncol = 365)
  simu_x<-t(t(simu_x)-colMeans(simu_x))
  saveRDS(simu_x, file = paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Rainfall_data_processed/grid",gridcode,"_rain.rds"))
  
  
  if(gridcode<=1600) Streamflow<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/StreamflowData/Lab-Wide/", gridcode,"_fabric_hydroCAN_2021data.csv"))
  if(gridcode>1600) Streamflow<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/StreamflowData/Lab-Wide/", gridcode,"_fabric_hydroUS_combinedData_revJune2021.csv"))
  Streamflow<-Streamflow[Streamflow$Year>=1979 & Streamflow$Year<2019,]
  
  Streamflow$date<-as.Date(paste0(Streamflow$Year,"-",Streamflow$Month,"-",Streamflow$Day))
  
  DateSeq<-seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day")
  #DateSeq2<-seq(from=as.Date("1979-10-01"),to=as.Date("2018-09-30"),by="day")
  
  dateDat<-data.frame(date=DateSeq)
  
  Fullstream<-merge(dateDat,Streamflow,by="date",all.x=T)
  Streamflow<-Fullstream$Streamflow_mm.d
  Streamflow[Streamflow<0]<-NA
  Streamflow<-log(Streamflow+1)
  leapdays<-which(month(DateSeq)==2 & day(DateSeq)==29)
  Streamflow[leapdays+1]<-(Streamflow[leapdays] + Streamflow[leapdays+1])/2
  Streamflow<-Streamflow[-leapdays]
  
  df<-data.frame(y=Streamflow,R=Rain,R1=c(NA,Rain[1:(length(Rain)-1)]),R2=c(NA,NA,Rain[1:(length(Rain)-2)]),S1=c(NA,Streamflow[1:(length(Streamflow)-1)]),
                 S2=c(NA,NA,Streamflow[1:(length(Streamflow)-2)]),doy=rep(1:365,length(Streamflow)/365))
  
  df<-missRanger(df)
  
  
  
  Streamflow<- df$y #fill in missing data
  
  Y=matrix(Streamflow,ncol = ncol(simu_x),byrow = T)
  Y<-t(t(Y)-colMeans(Y)) # subtract off alpha(t)
  saveRDS(Y, file = paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/streamflow_data_processed/grid",gridcode,"_streamflow.rds"))
}


rm(Daily_tmean,DailyTotalprcp)


###########################################################################################################
######################################   process data   ###################################################
############################################################################################################

#did 1344, 1405, 1409, 134, 303, 320, 355, 2052, 2042, 1795, 1759
#do 1344, 1759
# 1405, 320, and 2042 has differences in seasons
gridcode<-1759

simu_x<-readRDS(file=paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Rainfall_data_processed/grid",gridcode,"_rain.rds"))
Y<-readRDS(file = paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/streamflow_data_processed/grid",gridcode,"_streamflow.rds"))

M = ncol(simu_x)-1
max_lag<-120
K = (M+1)*max_lag
n_hp<-65

#regularization
RegMat<-get_D2(M,max_lag)
hProd<-crossprod(RegMat[[1]])
vProd<-crossprod(RegMat[[2]])

x<-as.numeric(t(simu_x))
y<-as.numeric(t(Y))
y<-y[max_lag:length(y)]

trainfrac<-0.6
valfrac<-0.2
testfrac<-0.2


y_train<-y[1:round(length(y)*trainfrac)]
y_val<-y[(round(length(y)*trainfrac)+1):round(length(y)*(trainfrac+valfrac))]
y_test<-y[(round(length(y)*(trainfrac+valfrac))+1):length(y)]
y_full<-y[1:round(length(y)*(trainfrac+valfrac))]

x_train<-x[1:(max_lag+length(y_train)-1)]
x_val<-x[(length(y_train)+1):(max_lag+length(y_train)+length(y_val)-1)]
x_test<-x[(length(y_train)+length(y_val)+1):length(x)]
x_full<-x[1:(max_lag+length(y_train)+length(y_val)-1)]

first_day<-1
first_day_train<-1
first_day_val<-(length(y_train)+1) %% 365
first_day_test<-(length(y_train)+length(y_val)+1) %% 365
first_day_full<-1

Xmat_train<-get_x(x_train,M,max_lag,first_day_train)
Xmat_val<-get_x(x_val,M,max_lag,first_day = first_day_val)
Xmat_test<-get_x(x_test,M,max_lag,first_day = first_day_test)
Xmat_full<-get_x(x_full,M,max_lag,first_day = first_day_full)
Xmat<-get_x(x,M,max_lag,first_day = first_day)

XXProd_train<-crossprod(Xmat_train)
XXProd_val<-crossprod(Xmat_val)
XXProd_test<-crossprod(Xmat_test)
XXProd_full<-crossprod(Xmat_full)
XXProd<-crossprod(Xmat)

XyProd_train<-crossprod(Xmat_train,matrix(y_train,ncol = 1))
XyProd_val<-crossprod(Xmat_val,matrix(y_val,ncol = 1))
XyProd_full<-crossprod(Xmat_full,matrix(y_full,ncol = 1))
XyProd_test<-crossprod(Xmat_test,matrix(y_test,ncol = 1))
XyProd<-crossprod(Xmat,matrix(y,ncol = 1))

############################################################################################################################################3
################################################    test on real streamflow data    #########################################################
#############################################################################################################################################

To_Zero<-rep(F,K)
result<-optimize_hp(hp_range=list(c(10,25),c(-12,15)),XTX=XXProd_train,hTh=hProd,vTv=vProd,XTy=XyProd_train,Xval=Xmat_val
                    ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)


best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
plot(result[[1]],result[[3]])
plot(result[[2]],result[[3]])

p_over_n<-sum(!To_Zero)/length(y_full)
#p_over_n<-length(y_train)/sum(!To_Zero)
#p_over_n<-1
pred_beta <- as.numeric(Matrix::solve(XXProd_full+ w_h[best]*p_over_n*hProd + w_v[best]*p_over_n*vProd, XyProd_full, sparse = TRUE))
plot_bvec_better(pred_beta,M,max_lag)

#detect where low betas are
grpnorms<-get_grp_normsH(pred_beta,M,max_lag)
plot_bvec_better(grpnorms,M,max_lag)


#optimize quantile for thresholding
q_result<-optimize_q(grpnorms = grpnorms,Nq=500,Xmat_train = Xmat_full,Xmat_val = Xmat_full,
                     y_train = y_full,y_val = y_full,hMat = RegMat[[1]],vMat = RegMat[[2]],w_h_best = w_h[best],w_v_best = w_v[best])
q_ngn<-q_result[[1]]
R2<-q_result[[2]]
plot(q_ngn,R2)
# distance2<-get_elbow_menger(q_ngn,R2)
# whichbestq_menger<-which.max(distance2[[3]])
# finalq<-distance2[[1]][whichbestq_menger]


distances1<-get_elbow(q_ngn,R2)
whichbestq<-which.max(distances1)
points(q_ngn[whichbestq],R2[whichbestq],col="red")
q_ngn2<- q_ngn[q_ngn<= q_ngn[whichbestq]]
R22<-R2[q_ngn<=q_ngn[whichbestq]]
distances2<-get_elbow(q_ngn2,R22)
plot(q_ngn2,R22)
whichbestq2<-which.max(distances2)
points(q_ngn2[whichbestq2],R22[whichbestq2],col="red")
finalq<-q_ngn2[whichbestq2]
if(max(distances2)<0.2) finalq<-q_ngn[whichbestq]

print(finalq)
To_Zero<-grpnorms<quantile(grpnorms,finalq)

#redo smoothness
XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])

result<-optimize_hp(hp_range=list(c(10,25),c(-12,15)),XTX=XXProd_sparse,hTh=hProd_sparse,vTv=vProd_sparse,XTy=XyProd_sparse,Xval=Xmat_val
                    ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)

best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
#refit model and update beta
XyProd_sparse<-crossprod(Xmat_full[,!To_Zero],matrix(y_full,ncol = 1))
XXProd_sparse<-crossprod(Xmat_full[,!To_Zero])
hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])

p_over_n<-sum(!To_Zero)/length(y_full)
#p_over_n<-length(y_full)/sum(!To_Zero)
p_over_n<-1
pred_beta<-rep(0,length(pred_beta))
pred_beta[!To_Zero]<-as.numeric(Matrix::solve(XXProd_sparse+ w_h[best]*p_over_n*hProd_sparse + w_v[best]*p_over_n*vProd_sparse, XyProd_sparse, sparse = TRUE))

Y_hat_test=Xmat_test %*% pred_beta
1-sum((y_test- Y_hat_test)^2)/sum((y_test- mean(y_full))^2)
plot_bvec_better(pred_beta,M,max_lag)



############################################################################################################################################3
################################################    final prediction    #########################################################
#############################################################################################################################################


set.seed(111)
To_Zero<-rep(F,K)
result<-optimize_hp(hp_range=list(c(10,25),c(-12,15)),XTX=XXProd_full,hTh=hProd,vTv=vProd,XTy=XyProd_full,Xval=Xmat_test
                    ,yval=y_test,ytrain = y_full,n_hp=n_hp,To_Zero = To_Zero)


best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
plot(result[[1]],result[[3]])
plot(result[[2]],result[[3]])

p_over_n<-sum(!To_Zero)/length(y)
#p_over_n<-1
pred_beta <- as.numeric(Matrix::solve(XXProd+ w_h[best]*p_over_n*hProd + w_v[best]*p_over_n*vProd, XyProd, sparse = TRUE))
plot_bvec_better(pred_beta,M,max_lag)

#detect where low betas are
grpnorms<-get_grp_normsH(pred_beta,M,max_lag)
plot_bvec_better(grpnorms,M,max_lag)


#optimize quantile for thresholding
q_result<-optimize_q(grpnorms = grpnorms,Nq=500,Xmat_train = Xmat,Xmat_val = Xmat,
                     y_train = y,y_val = y,hMat = RegMat[[1]],vMat = RegMat[[2]],w_h_best = w_h[best],w_v_best = w_v[best])
q_ngn<-q_result[[1]]
R2<-q_result[[2]]
plot(q_ngn,R2)

distances1<-get_elbow(q_ngn,R2)
whichbestq<-which.max(distances1)
points(q_ngn[whichbestq],R2[whichbestq],col="red")
q_ngn2<- q_ngn[q_ngn<= q_ngn[whichbestq]]
R22<-R2[q_ngn<=q_ngn[whichbestq]]
distances2<-get_elbow(q_ngn2,R22)
plot(q_ngn2,R22)
whichbestq2<-which.max(get_elbow(q_ngn2,R22))
points(q_ngn2[whichbestq2],R22[whichbestq2],col="red")
finalq<-q_ngn2[whichbestq2]
if(max(distances2)<0.2) finalq<-q_ngn[whichbestq]

# distance2<-get_elbow_menger(q_ngn,R2)
# whichbestq_menger<-which.max(distance2[[3]])
# points(distance2[[1]][whichbestq_menger],distance2[[2]][whichbestq_menger],col="red")
# finalq<-distance2[[1]][whichbestq_menger]
print(finalq)
To_Zero<-grpnorms<quantile(grpnorms,finalq)


#redo smoothness
XyProd_sparse<-crossprod(Xmat_full[,!To_Zero],matrix(y_full,ncol = 1))
XXProd_sparse<-crossprod(Xmat_full[,!To_Zero])
hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])

result<-optimize_hp(hp_range=list(c(10,25),c(-12,15)),XTX=XXProd_sparse,hTh=hProd_sparse,vTv=vProd_sparse,XTy=XyProd_sparse,Xval=Xmat_test
                    ,yval=y_test,ytrain = y_full,n_hp=n_hp,To_Zero = To_Zero)

best<-which.max(result[[3]])
w_h<-exp(result[[1]])
w_v<-exp(result[[2]])
#refit model and update beta
XyProd_sparse<-crossprod(Xmat[,!To_Zero],matrix(y,ncol = 1))
XXProd_sparse<-crossprod(Xmat[,!To_Zero])
hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])

p_over_n<-sum(!To_Zero)/length(y)
#p_over_n<-1
pred_beta<-rep(0,length(pred_beta))
pred_beta[!To_Zero]<-as.numeric(Matrix::solve(XXProd_sparse+ w_h[best]*p_over_n*hProd_sparse + w_v[best]*p_over_n*vProd_sparse, XyProd_sparse, sparse = TRUE))

Y_hat_test=Xmat %*% pred_beta
1-sum((y- Y_hat_test)^2)/sum((y- mean(y))^2)
plot(as.numeric(Y_hat_test),y- as.numeric(Y_hat_test))
plot_bvec_better(pred_beta,M,max_lag)

jpeg(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/groundtruth",gridcode,".jpeg"),width = 12, height = 12, units = 'in',res = 600)
plot_bvec_better(pred_beta,M,max_lag)
dev.off()

write.csv(data.frame(pred_beta=pred_beta),paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/groundtruth",gridcode,".csv"),row.names = F)




gridcode<-1759
pred_beta<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/groundtruth",gridcode,".csv"))
pred_beta<-pred_beta$pred_beta


############################################################################################################################################3
################################################    run simulation study       #########################################################
#############################################################################################################################################

gridcode<-1344
desiredR2<-0.8


simu_x<-readRDS(file=paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Rainfall_data_processed/grid",gridcode,"_rain.rds"))
x<-as.numeric(t(simu_x))

M = ncol(simu_x)-1
max_lag<-120
K = (M+1)*max_lag

#regularization
RegMat<-get_D2(M,max_lag)
hProd<-crossprod(RegMat[[1]])
vProd<-crossprod(RegMat[[2]])

trainfrac<-0.8
valfrac<-0.2

Y<-readRDS(file = paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/streamflow_data_processed/grid",gridcode,"_streamflow.rds"))
y<-as.numeric(t(Y))
y<-y[max_lag:length(y)]

y_train<-y[1:round(length(y)*trainfrac)]
y_val<-y[(round(length(y)*trainfrac)+1):length(y)]

x_train<-x[1:(max_lag+length(y_train)-1)]
x_val<-x[(length(y_train)+1):length(x)]

first_day<-1
first_day_train<-1
first_day_val<-(length(y_train)+1) %% 365

Xmat_train<-get_x(x_train,M,max_lag,first_day_train)
Xmat_val<-get_x(x_val,M,max_lag,first_day_val)
Xmat<-get_x(x,M,max_lag,first_day)

XXProd_train<-crossprod(Xmat_train)
XXProd_val<-crossprod(Xmat_val)
XXProd<-crossprod(Xmat)

bvec<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/groundtruth",gridcode,".csv"))
bvec<-bvec$pred_beta
true_delta<-compute_delta(bvec,M,max_lag)

n_hp<-65
nsim<-20
beta_R2<-rep(0,nsim)
delta_bias<-rep(0,nsim)
delta_cor<-rep(0,nsim)
for(i in 1:nsim){
  y=Xmat %*% bvec
  y<-as.numeric(y)
  noiseVar<-(var(y)-desiredR2*var(y))/desiredR2
  set.seed(i)
  y=y+rnorm(length(y),sd=sqrt(noiseVar))
  y_train<-y[1:round(length(y)*trainfrac)]
  y_val<-y[(round(length(y)*trainfrac)+1):length(y)]
  
  XyProd_train<-crossprod(Xmat_train,matrix(y_train,ncol = 1))
  XyProd_val<-crossprod(Xmat_val,matrix(y_val,ncol = 1))
  XyProd<-crossprod(Xmat,matrix(y,ncol = 1))
  
  
  To_Zero<-rep(F,K)
  result<-optimize_hp(hp_range=list(c(10,20),c(-12,15)),XTX=XXProd_train,hTh=hProd,vTv=vProd,XTy=XyProd_train,Xval=Xmat_val
                      ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)
  
  best<-which.max(result[[3]])
  bestR2<-max(result[[3]])
  w_h<-exp(result[[1]])
  w_v<-exp(result[[2]])
  p_over_n<-sum(!To_Zero)/length(y)
  #p_over_n<-1
  pred_beta <- as.numeric(Matrix::solve(XXProd+ w_h[best]*p_over_n*hProd + w_v[best]*p_over_n*vProd, XyProd, sparse = TRUE))
  #plot_bvec_better(pred_beta,M,max_lag)
  print(paste0("non sparese w_h and w_v are: ",round(w_h[best],3)," and ",round(w_v[best],3)))
  
  #detect where low betas are
  grpnorms<-get_grp_normsH(pred_beta,M,max_lag)
  #plot_bvec_better(grpnorms,M,max_lag)
  
  q_result<-optimize_q(grpnorms = grpnorms,Nq=400,Xmat_train = Xmat,Xmat_val = Xmat,
                       y_train = y,y_val = y,hMat = RegMat[[1]],vMat = RegMat[[2]],w_h_best = w_h[best],w_v_best = w_v[best])
  q_ngn<-q_result[[1]]
  R2<-q_result[[2]]
  plot(q_ngn,R2)
  # distance2<-get_elbow_menger(q_ngn,R2)
  # whichbestq_menger<-which.max(distance2[[3]])
  # points(distance2[[1]][whichbestq_menger],distance2[[2]][whichbestq_menger],col="red")
  # finalq<-distance2[[1]][whichbestq_menger]
  # print(finalq)
  # To_Zero<-grpnorms<quantile(grpnorms,finalq)
  
  
  distances1<-get_elbow(q_ngn,R2)
  whichbestq<-which.max(distances1)
  points(q_ngn[whichbestq],R2[whichbestq],col="red")


  q_ngn2<- q_ngn[q_ngn<= q_ngn[whichbestq]]
  R22<-R2[q_ngn<=q_ngn[whichbestq]]
  distances2<-get_elbow(q_ngn2,R22)
  plot(q_ngn2,R22)
  whichbestq2<-which.max(get_elbow(q_ngn2,R22))
  points(q_ngn2[whichbestq2],R22[whichbestq2],col="red")
  finalq<-q_ngn2[whichbestq2]
  if(max(distances2)<0.2) finalq<-q_ngn[whichbestq]
  print(finalq)
  To_Zero<-grpnorms<quantile(grpnorms,finalq)
  
  #redo smoothness
  XyProd_sparse<-crossprod(Xmat_train[,!To_Zero],matrix(y_train,ncol = 1))
  XXProd_sparse<-crossprod(Xmat_train[,!To_Zero])
  hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
  vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])
  
  result<-optimize_hp(hp_range=list(c(10,20),c(-12,15)),XTX=XXProd_sparse,hTh=hProd_sparse,vTv=vProd_sparse,XTy=XyProd_sparse,Xval=Xmat_val
                      ,yval=y_val,ytrain = y_train,n_hp=n_hp,To_Zero = To_Zero)
  
  best<-which.max(result[[3]])
  w_h<-exp(result[[1]])
  w_v<-exp(result[[2]])
  #refit model and update beta
  XyProd_sparse<-crossprod(Xmat[,!To_Zero],matrix(y,ncol = 1))
  XXProd_sparse<-crossprod(Xmat[,!To_Zero])
  hProd_sparse<-crossprod(RegMat[[1]][,!To_Zero])
  vProd_sparse<-crossprod(RegMat[[2]][,!To_Zero])
  
  p_over_n<-sum(!To_Zero)/length(y)
  #p_over_n<-1
  pred_beta<-rep(0,length(pred_beta))
  pred_beta[!To_Zero]<-as.numeric(Matrix::solve(XXProd_sparse+ w_h[best]*p_over_n*hProd_sparse + w_v[best]*p_over_n*vProd_sparse, XyProd_sparse, sparse = TRUE))
  #plot_bvec_better(pred_beta,M,max_lag)
  print(paste0("sparese w_h and w_v are: ",round(w_h[best],3)," and ",round(w_v[best],3)))
  
  beta_R2[i]<-compute_beta_R2(pred_beta,bvec)
  est_delta<-compute_delta(pred_beta,M,max_lag)
  delta_bias[i]<-mean(est_delta-true_delta)
  delta_cor[i]<-cor(est_delta,true_delta)
  
  write.csv(data.frame(bR2=beta_R2,db=delta_bias,dc=delta_cor),
            paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/SimStudy",gridcode,"_R2_",desiredR2,".csv"),row.names = F)
  print(i)
}


resultSim<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/SimStudy",1344,"_R2_",0.4,".csv"))
mean(resultSim$bR2)
sd(resultSim$bR2)

mean(resultSim$db)
sd(resultSim$db)

mean(resultSim$dc)
sd(resultSim$dc)


resultSim<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/SimStudy",1344,"_R2_",0.8,".csv"))

mean(resultSim$bR2)
sd(resultSim$bR2)

mean(resultSim$db)
sd(resultSim$db)

mean(resultSim$dc)
sd(resultSim$dc)

resultSim<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/SimStudy",1759,"_R2_",0.4,".csv"))

mean(resultSim$bR2)
sd(resultSim$bR2)

mean(resultSim$db)
sd(resultSim$db)

mean(resultSim$dc)
sd(resultSim$dc)

resultSim<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/SimStudy",1759,"_R2_",0.8,".csv"))

mean(resultSim$bR2)
sd(resultSim$bR2)

mean(resultSim$db)
sd(resultSim$db)

mean(resultSim$dc)
sd(resultSim$dc)



###########################################################################################################
######################################   plot temperature and rainfall and streamflow  ####################
############################################################################################################


gridcodes<-c(1344, 1759)

#get rainfall data
DailyTotalprcp <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/ClimateData/ClimateData/DailyTotalprcp.csv")
Daily_tmean <- read.csv("C:/Users/joeja/Desktop/researchMasters/Asad_paper/ClimateData/ClimateData/Daily_tmean.csv")

#gridcodes<-as.numeric(DailyTotalprcp[,1])


for(i in 1:length(gridcodes)){
  gridcode<-gridcodes[i]
  
  Precip<-as.numeric(DailyTotalprcp[DailyTotalprcp$Allresults==gridcode,2:ncol(DailyTotalprcp)])
  Temp<-as.numeric(Daily_tmean[Daily_tmean$Allresults==gridcode,2:ncol(Daily_tmean)])
  Rain<-Precip
  Rain[Temp<=0]<-0
  Rain[Rain<0.1]<-0
  
  #toadd<-min(Rain[Rain!=0],na.rm = T)/2
  #Rain[Rain==0]<-toadd
  #Rain<-log(Rain)
  
  #remove leap days from Rain
  #Date_seq<-seq(from=as.Date("1979-10-01"),to=as.Date("2018-09-30"),by="day")
  #Rain<-Rain[seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day") %in% Date_seq]
  Date_seq<-seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day")
  leapdays<-which(day(Date_seq)==29 & month(Date_seq)==2)
  before_lds<-leapdays-1
  Rain[before_lds]<-(Rain[before_lds]+Rain[leapdays])/2
  Rain<-Rain[-leapdays]
  #Rain[Rain<0.1]<-0 #ask ali about this
  
  #put into matrix form
  simu_x<-matrix(Rain,byrow = T,ncol = 365)
  AvgRain<-colMeans(simu_x)
  AvgTemp<-colMeans(matrix(Temp[-leapdays],byrow = T,ncol = 365))
  
  if(gridcode<=1600) Streamflow<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/StreamflowData/Lab-Wide/", gridcode,"_fabric_hydroCAN_2021data.csv"))
  if(gridcode>1600) Streamflow<-read.csv(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/StreamflowData/Lab-Wide/", gridcode,"_fabric_hydroUS_combinedData_revJune2021.csv"))
  Streamflow<-Streamflow[Streamflow$Year>=1979 & Streamflow$Year<2019,]
  
  Streamflow$date<-as.Date(paste0(Streamflow$Year,"-",Streamflow$Month,"-",Streamflow$Day))
  
  DateSeq<-seq(from=as.Date("1979-01-01"),to=as.Date("2018-12-31"),by="day")
  #DateSeq2<-seq(from=as.Date("1979-10-01"),to=as.Date("2018-09-30"),by="day")
  
  dateDat<-data.frame(date=DateSeq)
  
  Fullstream<-merge(dateDat,Streamflow,by="date",all.x=T)
  Streamflow<-Fullstream$Streamflow_mm.d
  Streamflow[Streamflow<0]<-NA
  leapdays<-which(month(DateSeq)==2 & day(DateSeq)==29)
  Streamflow[leapdays+1]<-(Streamflow[leapdays] + Streamflow[leapdays+1])/2
  Streamflow<-Streamflow[-leapdays]
  
  df<-data.frame(y=Streamflow,R=Rain,R1=c(NA,Rain[1:(length(Rain)-1)]),R2=c(NA,NA,Rain[1:(length(Rain)-2)]),S1=c(NA,Streamflow[1:(length(Streamflow)-1)]),
                 S2=c(NA,NA,Streamflow[1:(length(Streamflow)-2)]),doy=rep(1:365,length(Streamflow)/365))
  
  df<-missRanger(df)
  
  
  
  Streamflow<- df$y #fill in missing data
  
  Y=matrix(Streamflow,ncol = ncol(simu_x),byrow = T)
  AvgStreamflow<-colMeans(Y)
  
  jpeg(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/AvgRainStream",gridcode,".jpeg"),width = 12, height = 8, units = 'in',res = 600)
  plotmax<-max(AvgRain)*1.1
  par(mar=c(5, 4, 4, 6) + 0.1)
  plot(1:365, AvgRain, axes=FALSE, ylim=c(0,plotmax), xlab="", ylab="", 
       type="l",col="black",lwd=2)
  axis(2, ylim=c(0,plotmax),col="black",las=1)  ## las=1 makes horizontal labels
  mtext("Rainfall (mm)",side=2,line=2.5,cex=1.4)
  box()
  
  
  plotmax<-max(AvgStreamflow)*1.1
  par(new=TRUE)
  plot(1:365, AvgStreamflow,  xlab="", ylab="", ylim=c(0,plotmax), 
       axes=FALSE, type="l",lwd=2, col="red")
  mtext("Streamflow (mm)",side=4,col="red",line=4,cex=1.4) 
  axis(4, ylim=c(0,plotmax), col="red",col.axis="red",las=1)
  axis(1,pretty(range(1:365),10))
  mtext("Day of the year",side=1,col="black",line=2.5,cex=1.4)  
  legend("topleft",legend=c("Rainfall","Streamflow"),
         text.col=c("black","red"),pch=c(15,15),col=c("black","red"),cex=1.4)
  dev.off()
}


rm(Daily_tmean,DailyTotalprcp)




############################################################################################################################################3
################################################    play with code       #########################################################
#############################################################################################################################################



n_hp<-25

result<-optimize_hp(hp_range=c(0,50000000),XTX=XXProd_train,hTh=hProd,vTv=vProd,XTy=XyProd_train,Xval=Xmat_val
                    ,yval=y_val,n_hp=n_hp,To_Zero = rep(F,K))


best<-which.max(result[[3]])
w_h<-result[[1]]
w_v<-result[[2]]


pred_beta <- as.numeric(Matrix::solve(XXProd_test+ w_h[best]*hProd + w_v[best]*vProd, XyProd_test, sparse = TRUE))
plot_bvec_better(pred_beta,M,max_lag)

pred_beta <- as.numeric(Matrix::solve(XXProd_train+ w_h[best]*hProd + w_v[best]*vProd, XyProd_train, sparse = TRUE))
plot_bvec_better(pred_beta,M,max_lag)

Y_hat_val=Xmat_train %*% pred_beta
plot(y_train,as.numeric(y_train- Y_hat_val))
plot(as.numeric(Y_hat_val),as.numeric(y_train- Y_hat_val))
cor(as.numeric(Y_hat_val),as.numeric(y_train- Y_hat_val))

matplot(t(Y_val),type="l")
matplot(t(simu_x_val),type = "l")

pred_beta <- as.numeric(Matrix::solve(XXProd_val+ w_h[best]*hProd + w_v[best]*vProd, XyProd_val, sparse = TRUE))
plot_bvec_better(pred_beta,M,max_lag)

Y_hat_val=Xmat_val %*% pred_beta
plot(y_val,as.numeric(y_val- Y_hat_val))
plot(as.numeric(Y_hat_val),as.numeric(y_val- Y_hat_val))
cor(as.numeric(Y_hat_val),as.numeric(y_val- Y_hat_val))
hist(as.numeric(y_val- Y_hat_val))
mean(as.numeric(y_val- Y_hat_val))



library(glmnet)
set.seed(1234)
MC <- 50
l_CV <- numeric(MC)
ns   <- numeric(MC)
ps   <- numeric(MC)
y_var<- numeric(MC)
for(i in 1:MC){
  n       <- round(runif(1, 30, 1000))
  ns[i]   <- n
  p       <- round(n * runif(1, 1, 10))
  ps[i]   <- p
  X       <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y       <- X %*% rnorm(p) + rnorm(n)
  y_var[i] <- var(y)
  l_CV[i] <- cv.glmnet(scale(X), y, alpha = 0)$lambda.min
  cat("\n", i, ": n =", n, "p =", p, "l =", l_CV[i])
}
plot(l_CV ~ I(ps / ns))
points((ps/ns)[y_var<median(y_var)],l_CV[y_var<median(y_var)],col="red")

dat<-data.frame(x=sort(q_ngn),y=R2[order(q_ngn)])

datalengths<-seq(from=50,by=1,to=length(dat$x))
for(i in 1:length(datalengths)){
  mod<-lm(y~x,data=dat[1:datalengths[i],])
}

plot(dat$x,dat$y)
fm2DNase1 <- nls(y ~ alpha0+alpha1*(x-x1)+alpha2*(x-x1)*tanh((x-x1)/1e-5)+alpha3*(x-x2)*tanh((x-x2)/1e-5),
                 data = dat,
                 start = list(alpha0=max(dat$y),alpha1=-5,alpha2=-5,x1=0.8,alpha3=-5,x2=0.9))
summary(fm2DNase1)
newx<-seq(from=0,to=1,length.out=1000)
preds1<-predict(fm2DNase1,data.frame(x=newx))
preds2<-predict(fm2DNase1,data.frame(x=dat$x))

plot(dat$x,dat$y)
lines(newx,preds1,col="red")

plot(dat$x,dat$y-preds2)


DNase1 <- subset(DNase, Run == 1)
## using conditional linearity
fm2DNase1 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(xmid = 0, scal = 1),
                 algorithm = "plinear")
summary(fm2DNase1)

## without conditional linearity
fm3DNase1 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1))
summary(fm3DNase1)

## using Port's nl2sol algorithm
fm4DNase1 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1),
                 algorithm = "port")
summary(fm4DNase1)


#numerical proof that we should stop at onset of nonlinearity
n<-600
p<-300
prob_sig<-0.5
niter<-1
savedR2<-matrix(0,ncol = p,nrow = niter)
desiredR2<-0.8
for(it in 1:niter){
  X<-matrix(rnorm(n*p),nrow=n)
  beta<-runif(p,min=1.5,max=3)
  sig_flag<-sample(c(0,1),p,replace=T,prob = c(prob_sig,1- prob_sig))
  beta[sig_flag==0]<-0
  
  y<- X %*% matrix(beta,ncol=1)
  noiseVar<-(var(y)-desiredR2*var(y))/desiredR2
  y<- y + rnorm(length(y),sd=sqrt(noiseVar))
  
  for(i in 1:p){
    mod<-lm.fit(x=as.matrix(X[,order(beta)[p:i]],nrow=n),y=y)
    tempR2<-rep(0,1)
    for(r in 1:length(tempR2)){
      curSamp<-sample(1:length(y),length(y)*0.9999)
      tempR2[r]<- 1-sum((y[curSamp]-mod$fitted.values[curSamp])^2)/sum((y[curSamp]-mean(y[curSamp]))^2)
    }
    savedR2[it,i]<-median(tempR2)
  }
  print(it)
}

distance1<-get_elbow((1:p)/p,colMeans(savedR2))
distance2<-get_elbow_menger((1:p)/p,colMeans(savedR2))


plot((1:p)/p,colMeans(savedR2))
abline(v=((1:p)/p)[which.max(distance1)],col="blue",lwd=2)
abline(v=prob_sig,col="red",lwd=2)
abline(v=distance2[[1]][which.max(distance2[[3]])],col="black",lwd=2,lty=3)


jpeg(paste0("C:/Users/joeja/Desktop/researchMasters/Asad_paper/Methods_paper/numericalproof_",prob_sig,"_R2_",desiredR2,".jpeg"),width = 12, height = 12, units = 'in',res = 600)
plot((1:p)/p,colMeans(savedR2))
abline(v=prob_sig,col="red",lwd=2)
abline(v=((1:p)/p)[which.max(distance1)],col="blue",lwd=2)
abline(v=distance2[[1]][which.max(distance2[[2]])],col="black",lwd=2,lty=3)
dev.off()



