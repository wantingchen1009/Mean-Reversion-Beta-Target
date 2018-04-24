###First, read the file
factormodel_data<-read.csv('factormodel_data.csv',sep=',',header=TRUE);
colnames(factormodel_data) <-c("Date","FXE_return","EWJ_return","GLD_return","QQQ_return","SPY_return","SHV_return","DBA_return","USO_return","XBI_return","ILF_return","GAF_return","EPP_return","FXZ_return","rm_rf","SMB","HML","rf");
ticker<-colnames(factormodel_data);
CovFactors<-cov(factormodel_data[,15:17])

#factor model
rf<-matrix(factormodel_data[,18], ncol=40, nrow=60);
rm_rf<-matrix(factormodel_data[,15], ncol=40, nrow=60);
SMB<-matrix(factormodel_data[,16], ncol=40, nrow=60);
HML<-matrix(factormodel_data[,17], ncol=40, nrow=60);
return_data<-factormodel_data[1:60,2:14]
result<-function(k)
{
  famamodel<-function(i)
  {
    summary(lm(return_data[,i]-rf[,k]~rm_rf[,k]+SMB[,k]+HML[,k]))
    fama<-matrix(0,4,13);
    fama[,i]<-summary(lm(return_data[,i]-rf[,k]~rm_rf[,k]+SMB[,k]+HML[,k]))$coefficients[,1]
    return(fama[,i])
  }
return_data_fama<-lapply(1:13,function(i){return(famamodel(i))})
return_data_fama_matrix<-matrix(unlist(return_data_fama),ncol = 4, byrow = TRUE)
return_data_residual<-lapply(1:13,function(i){return(residuals(lm(return_data[,i]-rf[,k]~rm_rf[,k]+SMB[,k]+HML[,k])))})
return_data_residual_matrix<-matrix(unlist(return_data_residual),ncol =60, byrow = TRUE)
Dcov <- apply(return_data_residual_matrix,1,var)
RES <-diag(c(Dcov))  
covRtMatrix<-return_data_fama_matrix[,2:4]%*%CovFactors%*%t(return_data_fama_matrix[,2:4])+RES##

FXE_P<-lapply(1:13,function(i){return((return_data[,i]-rep(return_data_fama_matrix[i,1],60)-return_data_residual[[i]])/60+rep(return_data_fama_matrix[i,1],60))})
FXE_P_matrix <- matrix(unlist(FXE_P),ncol=13,byrow=FALSE)

mylist<-list(fama=return_data_fama_matrix, resi_dual=RES, expectation=FXE_P_matrix,covRt=covRtMatrix)

return(mylist)

}
parameter<-lapply(1:40,function(k){return(result(k))})
#parameter[[1]]$resi_dual


#all the data 
# beta of security Si,  vector 1*13
CAPMbeta<-c() 
CAPMbetaFunction<- function(i){
  CAPMbeta[i-1]=cov(factormodel_data[,i],(factormodel_data[,15]+factormodel_data[,18]))/((sd(factormodel_data[,15]+factormodel_data[,18]))^2)
}
CAPMbeta <- unlist(lapply(2:14,CAPMbetaFunction))
# ExpectedReturn list 2000 1*13

ExpectedReturn <- rbind(parameter[[1]]$expectation,parameter[[2]]$expectation,parameter[[3]]$expectation,parameter[[4]]$expectation,parameter[[5]]$expectation,
                        parameter[[6]]$expectation,parameter[[7]]$expectation,parameter[[8]]$expectation,parameter[[9]]$expectation,parameter[[10]]$expectation,
                        parameter[[11]]$expectation,parameter[[12]]$expectation,parameter[[13]]$expectation,parameter[[14]]$expectation,parameter[[15]]$expectation,
                        parameter[[16]]$expectation,parameter[[17]]$expectation,parameter[[18]]$expectation,parameter[[19]]$expectation,parameter[[20]]$expectation,
                        parameter[[21]]$expectation,parameter[[22]]$expectation,parameter[[23]]$expectation,parameter[[24]]$expectation,parameter[[25]]$expectation,
                        parameter[[26]]$expectation,parameter[[27]]$expectation,parameter[[28]]$expectation,parameter[[29]]$expectation,parameter[[30]]$expectation,
                        parameter[[31]]$expectation,parameter[[32]]$expectation,parameter[[33]]$expectation,parameter[[34]]$expectation,parameter[[35]]$expectation,
                        parameter[[36]]$expectation,parameter[[37]]$expectation,parameter[[38]]$expectation,parameter[[39]]$expectation,parameter[[40]]$expectation)

library(quadprog)
#optimization 
#the first portfolio,target_beta=0.5
Optimization1Function1 <- function(i){
w1<- c()# my portfolio
n=13
lemda=0.1
target_beta=0.5
wp1_backtest<<- rep(1/n, n)
Dmat <- diag(rep(lemda,n)) 
dvec <- t(ExpectedReturn[i,]) 
Amat <- rbind(t(CAPMbeta),rep(1,n),diag(rep(1,n)),diag(rep(-1,n)))
bvec <- c(target_beta,1,rep(-2,n),rep(-2,n))
w <-  solve.QP(Dmat,dvec,t(Amat),bvec,meq=2,factorized=FALSE)
wp1_backtest <<- w$solution
}
Optimization1Wp1 <- lapply(1:2400,Optimization1Function1)
Optimization1Wp1Matrix<-matrix(unlist(Optimization1Wp1),ncol = 13, byrow = TRUE)
#the first portfolio,target_beta=1
Optimization1Function2 <- function(i){
  w1<- c()# my portfolio
  n=13
  lemda=0.1
  target_beta=1
  wp1_backtest<<- rep(1/n, n)
  Dmat <- diag(rep(lemda,n)) 
  dvec <- t(ExpectedReturn[i,]) 
  Amat <- rbind(t(CAPMbeta),rep(1,n),diag(rep(1,n)),diag(rep(-1,n)))
  bvec <- c(target_beta,1,rep(-2,n),rep(-2,n))
  w <-  solve.QP(Dmat,dvec,t(Amat),bvec,meq=2,factorized=FALSE)
  wp1_backtest <<- w$solution
}
Optimization1Wp2 <- lapply(1:2400,Optimization1Function2)
Optimization1Wp2Matrix<-matrix(unlist(Optimization1Wp2),ncol = 13, byrow = TRUE)

#the first portfolio,target_beta=1.5
Optimization1Function3 <- function(i){
  w1<- c()# my portfolio
  n=13
  lemda=0.1
  target_beta=1.5
  wp1_backtest<<- rep(1/n, n)
  Dmat <- diag(rep(lemda,n)) 
  dvec <- t(ExpectedReturn[i,]) 
  Amat <- rbind(t(CAPMbeta),rep(1,n),diag(rep(1,n)),diag(rep(-1,n)))
  bvec <- c(target_beta,1,rep(-2,n),rep(-2,n))
  w <-  solve.QP(Dmat,dvec,t(Amat),bvec,meq=2,factorized=FALSE)
  wp1_backtest <<- w$solution
}
Optimization1Wp3 <- lapply(1:2400,Optimization1Function3)
Optimization1Wp3Matrix<-matrix(unlist(Optimization1Wp3),ncol = 13, byrow = TRUE)


#the second  portfolio
Optimization2Function <- function(i){
w2<- c()# my portfolio
n=13
lemda=0.1
target_return=0.15
wp2_backtest <<- rep(1/n, n)
Dmat <- parameter[[ceiling(i/60)]]$covRt+ diag(rep(lemda,n)) 
dvec <- 2*lemda*wp2_backtest 
Amat <- rbind(t(ExpectedReturn[i,]),rep(1,n),diag(rep(1,n)),diag(rep(-1,n)))
bvec <- c(target_return,1,rep(-2,n),rep(-2,n))
w <-  solve.QP(Dmat,dvec,t(Amat),bvec,meq=2,factorized=FALSE)
wp2_backtest <<- w$solution
}
Optimization2Wp <- lapply(1:2400,Optimization2Function)
Optimization2WpMatrix<-matrix(unlist(Optimization2Wp),ncol = 13, byrow = TRUE)

###plot the value
value1<-c()
value2<-c()
value3<-c()
value4<-c()
value5<-c()
DailyReturn <- as.matrix(factormodel_data[,2:14],ncol=13, nrow=2400)
value1[1]<-t(Optimization1Wp1Matrix[1,])%*%DailyReturn[1,]
value2[1]<-t(Optimization1Wp2Matrix[1,])%*%DailyReturn[1,]
value3[1]<-t(Optimization1Wp3Matrix[1,])%*%DailyReturn[1,]
value4[1]<-t(Optimization2WpMatrix[1,])%*%DailyReturn[1,]
value5[1]<-factormodel_data[1,6]
PlotDataFunction <- function(i){
  value1[i] <<- t(Optimization1Wp1Matrix[i,])%*%DailyReturn[i,]+value1[i-1]
  value2[i] <<- t(Optimization1Wp2Matrix[i,])%*%DailyReturn[i,]+value2[i-1]
  value3[i] <<- t(Optimization1Wp3Matrix[i,])%*%DailyReturn[i,]+value3[i-1]
  value4[i] <<- t(Optimization2WpMatrix[1,])%*%DailyReturn[i,]+value4[i-1]
  value5[i] <<- factormodel_data[i,6]+value5[i-1]
}
PlotData <-lapply(2:2400,PlotDataFunction)

x_Date<- seq(from=factormodel_data[,1])
PlotDateFrame<- data.frame(x_Date,value1,value2,value3,value4,value5)
colnames(PlotDateFrame) <-c("MyDate","beta0.5","beta1","beta1.5","MinVol","SPY")
ticker<-colnames(PlotDateFrame);

plot(x_Date,PlotDateFrame$beta0.5,xaxt="n",type="l",col="black",main="cumulated PnLs",xlab="date",ylab="return(%)", ylim=c(-50, 200))
axis(1,at=1:2400,labels=factormodel_data[,1])
lines(x_Date,PlotDateFrame$beta1,col="red")
lines(x_Date,PlotDateFrame$beta1.5,col="blue")
lines(x_Date,PlotDateFrame$MinVol,col="yellow")
lines(x_Date,PlotDateFrame$SPY,col="green")
legend("topleft",cex=0.6,c("beta0.5","beta1","beta1.5","MinVol","SPY"),
       lty=c(1,1,1,1,1),col=c("black", "red","blue","yellow","green"))


###define a function to caculate the performance
library(PerformanceAnalytics)
ETF.returns<-read.csv('factormodel_data.csv',sep=',',header=TRUE);
coredata <- ETF.returns[2:14]
returndate <- as.POSIXlt(ETF.returns[,1])
rownames(coredata) <- returndate
returndataxts <- as.xts(coredata)
colnames(returndataxts) <- c("FXE","EWJ","GLD","QQQ","SPY","SHV","DBA","USO","XBI","ILF","GAF","EPP","FEZ")
r=1/100*mean(ETF.returns[,18])
MD<-function(x){
  round(maxDrawdown(x),10)
}  
SP<-function(x){
  sqrt(250)*SharpeRatio(x, Rf=r, FUN="StdDev")
}
perfm <- function(x)
{
  other.factor <- rbind(apply(x,2,sum),apply(x,2,MD),apply(x,2,SP))
  rownames(other.factor) <- c("Cumulated PnL","Max 10 days Drawdown","Sharpe Ratio")
  stat<-table.Stats(x)
  other.stat<-rbind(stat[7,],stat[3,],sqrt(250)*stat[14,],stat[15:16,])
  risk<-table.DownsideRisk(x)
  other.risk<-rbind(sqrt(250)*risk[10:11,])
  f<-rbind(other.factor,other.stat,other.risk)
  f
}


#######caculate the portfolio performance#######

library(PerformanceAnalytics)
#transfer the return to time series data
ETF.returns<-read.csv('factormodel_data.csv',sep=',',header=TRUE);
coredata <- 1/100*ETF.returns[2:14]
returndate <- as.POSIXlt(ETF.returns[,1])
rownames(coredata) <- returndate
returndataxts <- as.xts(coredata)
colnames(returndataxts) <- c("FXE","EWJ","GLD","QQQ","SPY","SHV","DBA","USO","XBI","ILF","GAF","EPP","FEZ")
rf=1/100*mean(ETF.returns[,18])

#define a function to caculate performance(annualized)
MD<-function(x){
  round(maxDrawdown(x),10)
}  
SP<-function(x){
  sqrt(250)*SharpeRatio(x, Rf=rf, FUN="StdDev")
}
perfm <- function(x)
{
  other.factor <- rbind(apply(x,2,sum),apply(x,2,MD),apply(x,2,SP))
  rownames(other.factor) <- c("Cumulated PnL","Max 10 days Drawdown","Sharpe Ratio")
  stat<-table.Stats(x)
  other.stat<-rbind(stat[7,],stat[3,],sqrt(250)*stat[14,],stat[15:16,])
  risk<-table.DownsideRisk(x)
  other.risk<-rbind(sqrt(250)*risk[10:11,])
  f<-rbind(other.factor,other.stat,other.risk)
  f
}
SPY<-returndataxts[,5]
colnames(SPY)<-NULL
spy<-as.numeric(SPY)

#apply the weight to ETF to structure the portfolio every day
mat1<-t(Optimization1Wp1Matrix)
mat2<-t(Optimization1Wp2Matrix)
mat3<-t(Optimization1Wp3Matrix)
mat4<-t(Optimization2WpMatrix)
rownames(coredata)<-NULL
colnames(coredata)<-NULL
dailyR<-as.matrix(coredata)
product1<-dailyR%*%mat1
product2<-dailyR%*%mat2
product3<-dailyR%*%mat3
product4<-dailyR%*%mat4
diag1<-1/100*diag(product1)
diag2<-1/100*diag(product2)
diag3<-1/100*diag(product3)
diag4<-1/100*diag(product4)
PortfolioReturn<-as.data.frame(cbind(diag1,diag2,diag3,diag4,spy))
rownames(PortfolioReturn) <- returndate
Portfolioxts <- as.xts(PortfolioReturn)
colnames(Portfolioxts) <- c("Target beta=0.5","Target beta=1","Target beta=1.5","Target return=15%","S&P 500")

##show the performance of 13 ETF and 4 portfolio##
perfm(returndataxts)
#portfolio before the financial crisis 3/26/2007 ~ 8/8/2007
perfm(Portfolioxts[1:95,])
#portfolio during the financial crisis 8/9/2007 ~ 7/1/2010
perfm(Portfolioxts[96:825,])
#portfolio after the financial crisis 7/2/2010 ~ 10/3/2016
perfm(Portfolioxts[826:2400,])

