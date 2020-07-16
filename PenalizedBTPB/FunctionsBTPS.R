#####Functions for Bivariate Tensor Product B-splines

##Required Libraries

require(mgcv)
require(msm)
require(stats)
require(rgl)
require(MASS)
require(reshape2)
require(plyr)
require(ggplot2)
########################Function Names and purpose
#fitBTPS-Fits a penalized bivariate tensor product B-spline using According to Section 3
#         Then estimates the first and second partial derivative with respect to the x values
#
#yieldDataSim-creates simulated data in accordance to section 4
#
#compute_bandwidth-Calculate bandwidth for kernel density estimation
#                   Helper function for VarianceEstimator
#
#predictDensity-calculate a bivariate kernel density, helper function for VarianceEstimator, uses compute_bandwidth
#
#VarianceEstimator-Calculates the variance for a fitBTPS object according to equation 22
#
#kernelDensity
#
#fitKernel

#######################################



##inputs: y- the dependent variable for the penalized BTPS
#x-the independent variable that will have it's derivative taken
#z-the second independent variable
#knots-the number of knots used for each axis
#degree-the degree of the spline used
#penalty the order of the penalty used
#tol1st- the tolerance when taking the first derivative
#tol2nd-the tolerance when taking the second derivative
#degree.sd-degree of the spline used when estimating sigma
#penalty.sd=c(1,1)-the penalty used when estimating sigma
#newXLim-the range of x values used to predict on
#newZLim-the range of z values used to predict on

##Returns a list with the following:
##outputs-origY a vector of the original y values used,
#origX-a vector of the original x values used
#origZ-a vector of the original z values used
#fitOrig-fitted values for the original data points
#newX-new x values to predict
#newZ-new z values to predict
#Fit-the predicted y of the new x and z values
#firstDerv- the estimated first derivatives at the new x and z values
#secondDerv- the estimated second derivatives at the new x and z values
#int.knots-the number of interior knots used for each axis
#smoothedSigma-the estimated sigma value, used for estimating the variance
#model-the output from the gam function

fitBTPS<-function(y,x,z,knots=c(0,0),degree=c(5,3),penalty=c(1,1),tol1st=.0000001,tol2nd=.1,degree.sd=c(3,3),penalty.sd=c(1,1),newXLim=c(min(x),max(x)),newZLim=c(min(z),max(z))){
  
  
  ##Get a value for the number of B-spline Basis needed for the set up based on degree and knots
  myK<-c(knots+degree+2)
  
  ##create the degree and penalty list for the function
  m<-list(c(degree[1],penalty[1]),c(degree[2],penalty[2]))

  ##fit the BTPS using the gam function in r
  theFit<-gam(y~te(x,z,bs=c("ps","ps"),m=m,k=myK),drop.intercept=F)
  
  ###create new data and fit the new data
  xNew<-seq(newXLim[1],newXLim[2],length.out = 50)##values of length.out might not work on all computers due to memory limitations
  zNew<-seq(newZLim[1],newZLim[2],length.out = 50)
  

  
  hold.dataframe<-expand.grid(xNew,zNew)
  
  xNew<-hold.dataframe$Var1
  zNew<-hold.dataframe$Var2

  
  yPred<-as.numeric(predict(theFit,data.frame(x=x,z=z)))
  Fit<-as.numeric(predict(theFit,data.frame(x=xNew,z=zNew)))
  
  ##Estimate the First derivative fit
  Derv1Fit<--(predict(theFit,data.frame(x=xNew,z=zNew))-predict(theFit,data.frame(x=xNew+tol1st,z=zNew)))/tol1st
  
  ##Estimate the Second derivative fit
  Derv2Fit<-(predict(theFit,data.frame(x=xNew-tol2nd,z=zNew))-2*predict(theFit,data.frame(x=xNew,z=zNew))+predict(theFit,data.frame(x=xNew+tol2nd,z=zNew)))/tol2nd^2
  
  
  ###Now for the error estimation
  
  Residuals.V<-as.numeric(theFit$residuals)^2
  m.sd<-list(penalty.sd,degree.sd)
  
  ##Fit the residuals as in equation 20
  fittedRes<-gam(Residuals.V~te(x,z,bs=c("ps","ps"),m=m.sd,k=myK),drop.intercept=F)
  residualsEst<-abs(as.numeric(predict(fittedRes,data.frame(x=xNew,z=zNew))))
  lambdaValues<-theFit$sp
  
  
  
  
  ##data returned
  return(list(origY=y,origX=x,origZ=z,fitOrig=yPred,newX=xNew,newZ=zNew,Fit=Fit,firstDerv=as.numeric(Derv1Fit),secondDerv=as.numeric(Derv2Fit),int.knots=knots,smoothedSigma=residualsEst,model=theFit))
  
}

##inputs: n- the number of data points to simulate
#YieldMeanFunc-the mean structure to use for the simulation
      #"Linear": (-25+1.3*z) "Quad":  1.2*(z-50)+(z-150)^2/200
#YieldVarFunc-the variance of the mean structure distribution
      #"Const":1 "NonConst":|z|^.2
#YieldError-the distribution of the yield structure
    #"Normal":10N(0,1) or "Beta":50[Beta(alpha,Beta)-alpha/(alpha+beta)]
#alphaE-alpha value if the beta yield error is used
#betaE-beta value if the beta yield error is used
#p-value used for p in equation 2
#meas_error-the measurement error added to the estimated premium

##Returns a list with the following:
##Yield_Curr-simulated yield values
##Yield_Hist-simulated historical yield values (APH)
##Cov_Rate-simulated coverage rate values
##Premium-simulated premium values
##Obs_Prem-premium values with the simulated measurement error
##OverallSigma-sigma values for the simulation


yieldDataSim<-function(n,YieldMeanFunc="Linear",YieldVarFunc="Const",YieldError="Normal",alphaE=5,betaE=3,p=1,meas_error=.001){
  
  ####simulate x and z values
  xSim<-sample(seq(.55,.95,by=.05),n,replace=T)
  zSim<-rtnorm(n,200,25,100,300)

  ##create the yield mean based on linear or quadratic
  if(YieldMeanFunc=="Linear"){
    yieldMean<--25+1.3*zSim
  }
  
  if(YieldMeanFunc=="Quad"){
    
    yieldMean<-1.2*(zSim-50)+(zSim-150)^2/200
  }
  
  ##create the sigma based on non-constant and constant
  if(YieldVarFunc=="NonConst"){
    yieldSigma<-abs(yieldMean)^.2
    
  }
  
  if(YieldVarFunc=="Const"){
    
    yieldSigma<-rep(1,n)
  }
  
  ###input the YieldError (Normal or Beta) and estimate the premium value 
  premium<-rep(NA,n)
  if(YieldError=="Normal"){
    error<-rnorm(n,0,1)
    
    
    if(YieldMeanFunc=="Linear"){s.factor<-10}
    if(YieldMeanFunc=="Quad"){s.factor<-10}
  
    myIntFunc<-function(x,mu,sigma){
      
      
      return(pnorm((x-mu)/sigma,0,1))
      
    }
    
    
    for(i in 1:n){
      premium[i]<-p*integrate(myIntFunc,lower=zSim[i]+yieldSigma[i]*s.factor*-5,upper=xSim[i]*zSim[i],mu=yieldMean[i],sigma=yieldSigma[i]*s.factor)$value
    }
    
  }
  
  if(YieldError=="Beta"){
    nonCentral<-alphaE/(alphaE+betaE)
    error<-rbeta(n,alphaE,betaE)-nonCentral

    if(YieldMeanFunc=="Quad"){s.factor<-50}
    if(YieldMeanFunc=="Linear"){s.factor<-50}
    #print(s.factor)
    
    myIntFunc<-function(x,alpha,beta,mu,sigma){
      
      
      return(pbeta(alpha/(alpha+beta)+(x-mu)/sigma,alpha,beta))
      
    }
    bottomInt<--alphaE/(alphaE+betaE)
    
    for(i in 1:n){
      premium[i]<-p*integrate(myIntFunc,lower=zSim[i]+yieldSigma[i]*s.factor*bottomInt,upper=xSim[i]*zSim[i],mu=yieldMean[i],sigma=yieldSigma[i]*s.factor,alpha=alphaE,beta=betaE)$value
      
      
      
    }
  }
  
  ##Calculate the current yield
  
  currentYield<-yieldMean+yieldSigma*s.factor*error


  ###Return List
  return(list(Yield_Curr=currentYield,Yield_Hist=zSim,Cov_Rate=xSim,Premium=premium,Obs_Prem=premium+rnorm(length(premium),0,meas_error),OverallSigma=yieldSigma*s.factor))
}

#####

##inputs: bw- a charter string of the bandwidth selector to use for gaussian kernels
#x-the values of which to estimate the bandwidth for


## outputs:returns vector of bandwidths evaluated at x
compute_bandwidth <- function(bw, x) {
  if (is.numeric(bw)) return(bw)
  if (is.character(bw)) bw <- match.fun(bw)
  
  if (is.function(bw)) {
    bw <- bw(x)
    message("Using bandwidth ", format(bw, digits = 3))
    bw
  } else {
    stop("Invalid bandwidth")
  }
}

##inputs: xOrig-x values used for fitting a bivariate kernel density
#yOrig-y values used for fitting a bivariate kernel density
#(xOrig, yOrig) are a paired observation
#xNew-x values on which to evaluate the fitted a bivariate kernel density
#yNew-y values on which to evaluate the fitted a bivariate kernel density
#(xNew, yNew) are a paired observation

## outputs: returns a list with a vector of density values for the (xNew, yNew) observations
predictDensity <- function(xOrig, yOrig,xNew,yNew) {
  
  n_out <- length(xNew)
  n_in<-length(xOrig)
  xh <- suppressMessages(compute_bandwidth("bw.nrd", xOrig))
  yh <- suppressMessages(compute_bandwidth("bw.nrd", yOrig))
  ax <- outer(xNew, xOrig, "-") / xh
  ay <- outer(yNew, yOrig, "-") / yh
  myDensity<- rowSums(matrix(stats::dnorm(ax) * stats::dnorm(ay), n_out, n_in) * 
                        1 / (length(xOrig) * xh * yh))
  
  return(list(myDensity))
}

##inputs: BTPS.obj-x values used for fitting a bivariate kernel density
#Derv: in vector form, the partial derivatives the variance should be calculated for 

## outputs:a list with the following values
#xValues: a vector of the corresponding x values of the variance values
#zValues: a vector of the corresponding z values of the variance values
#variance: a vector of the estimated variance for the points (xValues,zValues)

VarianceEstimator<-function(BTPS.obj,Derv=c(2,0),CropData=T){

  #calculate n_1 and n_1 as well as knots
  n_ind<-c(length(unique(BTPS.obj$origX)),length(BTPS.obj$origY)/length(unique(BTPS.obj$origX)))
  knots<-BTPS.obj$int.knots
  intHvalues<-c(NA,NA)
  hValues<-c(NA,NA)
  
  #get the needed info from the gam object
  penalties<-c( BTPS.obj$model$smooth[1][[1]]$margin[[1]]$null.space.dim, BTPS.obj$model$smooth[1][[1]]$margin[[2]]$null.space.dim)
  smoothingPenalties<-c(max(as.numeric(BTPS.obj$model$sp[1]),1),max(as.numeric(BTPS.obj$model$sp[2]),1))
  for(i in 1:2){
    
    
   # function for calculating H as seen in appendix A
    Hfunction<-function(x,dervs=2,penalty=2){
      if(penalty==1){equ<-expression((1/(2)*exp(-x)), "x")}
      if(penalty==2){equ<-expression(1/(2*sqrt(2))*exp(-1/(sqrt(2))*x)*(sin(x/sqrt(2))+cos(x/sqrt(2))), "x")}
      if(penalty==3){equ<-expression(1/(6)*exp(-x)+exp(-1/(2)*x)*(cos(x*sqrt(3)/2)/6+sin(x*sqrt(3)/2)*sqrt(3)/6), "x")}
      if(penalty==4){equ<-expression(exp(-.9239*x)*(.231*cos(.3827*x)+.0957*sin(.3827*x))+exp(-.3827*x)*(.0957*cos(.9239*x)+.231*sin(.9239*x)), "x")}
      
      if(dervs==0){equ<-equ[1]}
      if(dervs>0){
        for(i in 1:dervs){
          equ<-D(equ,"x")
          
        }
      }
      
      
      return(eval(equ)^2)
      
    }
    
    ##integration of the H function
    intHvalues[i]<-integrate(Hfunction,lower=.00001,upper=Inf,dervs=Derv[i],penalty=penalties[i])$value
    
    ##calculation of the h values as given in lemma 1
    hValues[i]<-((((smoothingPenalties[i]*knots[i]/n_ind[i]))^(1/(2*penalties[i])))/knots[i])^(2*Derv[i]+1)
  n<-1
    
    }
  finalH<-1

  #combine the h_1 and h_2 values into h
  if(CropData){
   finalH<-hValues[1]*hValues[2]

    
  }

  #estimate f hat by a bivariate kernal density estimator
  kernelDensity<-predictDensity(BTPS.obj$origX,BTPS.obj$origZ,BTPS.obj$newX,BTPS.obj$newZ)
  
  #combine into the final variance
  EstVar<-4*intHvalues[1]*intHvalues[2]*BTPS.obj$smoothedSigma/unlist(kernelDensity)/(finalH*length(BTPS.obj$origX))
  return(list(xValues=BTPS.obj$newX,zValues=BTPS.obj$newZ,variance=EstVar))
}




##inputs: w: the dependent variable on which to fit the spline
#z: the independent variable on which to fit the spline
#m: a vector of c(degree,penalty) used to fit the p-spline
#knots: number of interior knots to use 
#fixedZ:the z values at which to predict the fitted kernel estimator with
#predictW: the w values at which to predict the fitted kernel estimator with

## outputs: fHat a dataframe  of the fitted kernel estimator values at points  with columns fixedZ and predictW rows)

kernelEstimator<-function(w,z,m,knots,fixedZ,predictW){
  
  ##fit the p-spline to the data
  theFit<-gam(w~s(z,bs="ps",m=m,k=knots),drop.intercept=F)
  
  
  #fit the squared residuals to get sigma^2
  sigmaSquare<-gam(residuals(theFit)^2~s(z,bs="ps",m=m,k=knots),drop.intercept=F)
  
  
  #Get the estimated mean and sigma^2
  Predicted_Mu<-as.numeric(predict(theFit,data.frame(z=fixedZ)))
  
  Predicted_SigmaSquare<-abs(as.numeric(predict(sigmaSquare,data.frame(z=fixedZ))))
  
  
  
  
  
  
  epsilon<-residuals(theFit)/sqrt(predict(sigmaSquare))
  
  #Iterate through for each value of w and z
  fHat<-data.frame(predictW)
  names(fHat)<-"w"
  for(i in 1:length(fixedZ)){
    
    minVal<-(min(fHat$w,na.rm=T)-Predicted_Mu[i])/sqrt(Predicted_SigmaSquare[i])
    maxVal<-(max(fHat$w,na.rm=T)-Predicted_Mu[i])/sqrt(Predicted_SigmaSquare[i])
    ifBad<-tryCatch(density(epsilon,n=length(predictW),from=minVal, to=maxVal,na.rm=T),error=function(e){"ERROR"})
    if(ifBad[1]=="ERROR"){browser()}
    
    mydensity<-density(epsilon,n=length(predictW),from=minVal, to=maxVal,na.rm=T)
    mydensity$y
    
    fHat<-cbind(fHat,mydensity$y/sqrt(Predicted_SigmaSquare[i]))
    names(fHat)[ncol(fHat)]<-fixedZ[i]
    
  }
  
  return(fHat)
}

##inputs: w: the dependent variable on which to fit the spline
#z: the independent variable on which to fit the spline
#degree: the degree of penality used in the p-spline
#knots: number of interior knots to use 
#penalty: The penalty of the p-spline used
#fixedZ: the z values at which to predict the fitted kernel estimator with
#NumW: the number of w values to estimate the kernel density estimate at
#wRange: the range of w values to estimate the kernel density estimate at


## outputs: list with the following values
#estimate: The kernel density estimate
#SD: the standard error of the kernel density estimate
#fixedZvalue: THe z values that the kernel density was estimated at
fitKernel<-function(w,z,knots=1,degree=2,penalty=2,fixedZ=NA,NumW=1000,wRange=NA){
 
  #put knots kernel and degree and penalty into form to use gam function
   myK<-c(knots+degree+1)
  
  m<-c(degree,penalty)
  
  
  
  ####Kernel density estimation
  #set up the values at which to estimate at (default 5 z values)
  if(is.na(wRange[1])){
    
    new_w<-seq(min(w),max(w),length.out=NumW)
  }else{   new_w<-seq(min(wRange),max(wRange),length.out=NumW)}
  
  if(is.na(fixedZ[1])){
    
    predictionValues<-seq(min(z),max(z),length.out = 5)
    
    
  }else{predictionValues<-fixedZ}
  
  ##get the kernel density estimates
  finalFhat<-as.matrix(kernelEstimator(w,z,m,myK,predictionValues,new_w))

  
  ###Now for the error estimation using the jackknife method
  standardErrorMat<-array(,dim=c(NumW,length(predictionValues)+1,length(w)))
  for(i in 1:length(w)){

    standardErrorMat[,,i]<-(as.matrix(kernelEstimator(w[-i],z[-i],m,myK,predictionValues,new_w))-finalFhat)^2
    
  
  }
  
  standardErrorF<-apply(standardErrorMat, c(1,2), sum,na.rm=T)*(length(w)-1)/length(w)
  
  
  ##return the values
  return(list(estimate=finalFhat,SD=sqrt(standardErrorF),fixedZvalue=predictionValues))
  
}

