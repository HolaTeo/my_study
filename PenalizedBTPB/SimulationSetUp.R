##inputs: reps- the number of monte carlo replacations to run
#seed-the value to set the initial seed to
#n- the number of data points to simulate
#YieldMeanFunc-the mean structure to use for the simulation
#"Linear": (-25+1.3*z) "Quad":  1.2*(z-50)+(z-150)^2/200
#YieldVarFunc-the variance of the mean structure distribution
#"Const":1 "NonConst":|z|^.2
#YieldError-the distribution of the yield structure
#"Normal":10N(0,1) or "Beta":50[Beta(alpha,Beta)-alpha/(alpha+beta)]
#alphaE-alpha value if the beta yield error is used
#betaE-beta value if the beta yield error is used
#pVal-value used for p in equation 2
#meas_error-the measurement error added to the estimated premium
#BTPBknots-the number of knots used for each axis
#BTPBdegree-the degree of the spline used
#BTPBpenalty the order of the penalty used
#tol2nd-the tolerance when taking the second derivative
#newXLim-the range of x values used to predict on
#newZLim-the range of z values used to predict on

##Returns a list with the following:
#Estimates-a three dimensional matrix whose x axis is the coverage rate, y axis is the aph, and z is the replicate fit
#the values are the estimates for the given replicate, the first matrix in the rep dimension are the true values

#StdError-a three dimensional matrix with the standard error for each estimate, no true values are included


##Now the sim Study
simStudyFunc<-function(reps=1000,seed=9389,samplesize=1500,Yield_Mean="Linear",VarFunc="Const",Yield_Error="Normal"  ,alphaVal=5,betaVal=3,pVal=4,MeasurementError=.001,BTPBknots=c(2,9),BTPBpenalty=c(3,3),BTPBdegree=c(5,3),tol2nd = .00005,XLim =c(.55,.95) ,ZLim =c(100,300)){
  
  #set the seed
  set.seed(seed)
  
  #Create a list of the names
  repNames<-list("True")
 
  for(j in 2:(reps+1)){
    repNames[j]<-paste("Rep",(j-1),sep="")
    
  }
  
 
  Estimates<-array(,dim=c(50,50,(reps+1)))
  VarEst<-array(,dim=c(50,50,reps))
  

  ##iterate through the replicates
  for(reps in 2:(reps+1)){

    #sim the data
    sim.dat<-yieldDataSim(samplesize,YieldMeanFunc=Yield_Mean,YieldVarFunc=VarFunc,YieldError=Yield_Error,alphaE=alphaVal,betaE=betaVal,p=pVal,meas_error=MeasurementError)
    #fit the data
    test<-fitBTPS(sim.dat$Obs_Prem,sim.dat$Cov_Rate,sim.dat$Yield_Hist,knots=BTPBknots,penalty=BTPBpenalty,degree=BTPBdegree,tol2nd = tol2nd,newXLim =XLim ,newZLim =ZLim)
    
    mine<-VarianceEstimator(test) 
    overallSig<-rep(sim.dat$OverallSigma[1],length(test$newX))
    
    ##record the estimates
    Estimates[,,reps]<-matrix(test$secondDerv/(pVal*test$newZ^2),nrow=50)
    VarEst[,,reps-1]<-matrix(sqrt(mine$variance)/(pVal*test$newZ^2),nrow=50) 
    
    
    
    print(reps)
    
    
  }
  dimnames(Estimates)<-list(unique(test$newX),unique(test$newZ),repNames)
  dimnames(VarEst)<-list(unique(test$newX),unique(test$newZ),repNames[-1])
  
  #calculate the True Values (Depends on Set Up)
  {
    
    if(Yield_Mean=="Linear"){
      reg_Yield<-(-25+1.3*test$newZ)

    }
    
    if(Yield_Mean=="Quad"){
      reg_Yield<-(1.2*(test$newZ-50)+(test$newZ-150)^2/200)
    }
    
    if(Yield_Error=="Normal"){
    if(VarFunc=="Const"){
      overallSig<-rep(10,length(test$newX))
      
    }
      if(VarFunc=="NonConst"){
        overallSig<-10*abs(reg_Yield)^.2
        
      }  
      Estimates[,,1]<-matrix(dnorm((test$newX*test$newZ-reg_Yield)/overallSig)/overallSig,nrow=100) 
    }
    
    if(Yield_Error=="Beta"){
      if(VarFunc=="Const"){
        overallSig<-rep(50,length(test$newX))
        
      }
      if(VarFunc=="NonConst"){
        overallSig<-50*abs(reg_Yield)^.2
        
      }  
      Estimates[,,1]<-matrix(dbeta((test$newX*test$newZ-reg_Yield)/overallSig+(alphaVal/(alphaVal+betaVal)),alphaVal,betaVal)/overallSig,nrow=100)
      
    }
    
    

  

  
  
}

  
  return(list(Estimates=Estimates,Variance=VarEst))
}

##inputs: Nreps- the number of monte carlo replacations to run
#seed-the value to set the initial seed to
#samplesize- the number of data points to simulate
#zValues: the z values at which to predict the fitted kernel estimator with
#wRange: the range of w values to estimate the kernel density estimate at
#numberOfW: number of w values to estimate the kernel density at
#YieldMeanFunc-the mean structure to use for the simulation
#"Linear": (-25+1.3*z) "Quad":  1.2*(z-50)+(z-150)^2/200
#YieldVarFunc-the variance of the mean structure distribution
#"Const":1 "NonConst":|z|^.2
#YieldError-the distribution of the yield structure
#"Normal":10N(0,1) or "Beta":50[Beta(alpha,Beta)-alpha/(alpha+beta)]
#alphaE-alpha value if the beta yield error is used
#betaE-beta value if the beta yield error is used

##Returns a list with the following:
#Estimates-a three dimensional matrix whose x axis is the coverage rate, y axis is the aph, and z is the replicate fit
#the values are the estimates for the given replicate, the first matrix in the rep dimension are the true values

#Variance-a three dimensional matrix with the variance for each estimate, no true values are included



#Now the sim study for Kernel
simStudyKernel<-function(Nreps=1000,seed=9389,samplesize=1500,zValues=seq(100,300, by=10),wRange=c(0,450),numberOfW=1000,Yield_Mean="Linear",VarFunc="Const",Yield_Error="Normal"  ,alphaVal=5,betaVal=3){
  set.seed(seed)
  
  #set the data frames up to hold the simulated data and mc estimates
  repNames<-c()

  for(j in 1:Nreps){
    repNames[j]<-paste("Rep",(j),sep="")
    
  }
  fixedZvalue<-zValues
  numberFixed<-length(zValues)
  beginW<-wRange[1]
  endW<-wRange[2]
  Estimates<-array(,dim=c(numberOfW,numberFixed,Nreps+1),dimnames=list(seq(beginW,endW,length.out=numberOfW),fixedZvalue,c("True",repNames)))
  VarEst<-array(,dim=c(numberOfW,numberFixed,Nreps),dimnames=list(seq(beginW,endW,length.out=numberOfW),fixedZvalue,repNames))

  #interate through the replicates
  for(reps in 1:Nreps){
    print(paste("Replicate:", reps, " of ",Nreps))
 
    #sim the data
    sim.dat<-yieldDataSim(samplesize,YieldMeanFunc=Yield_Mean,YieldError=Yield_Error,YieldVarFunc = VarFunc)

    #fit the data
    fit<-fitKernel(sim.dat$Yield_Curr,sim.dat$Yield_Hist,knots=10,wRange=wRange,fixedZ = fixedZvalue,NumW = numberOfW)
    
#save the estimates
    Estimates[,,reps+1]<-fit[[1]][,-1]
    VarEst[,,reps]<-fit[[2]][,-1]
 

  }
  
##calculate the true values
  for(i in 1:numberOfW){
    thisW<-as.numeric(dimnames(Estimates)[[1]][i])
  
  #True Values (Depends on Set Up)
  {
    
    if(Yield_Mean=="Linear"){
      reg_Yield<-(-25+1.3*fixedZvalue)
      
    }
    
    if(Yield_Mean=="Quad"){
      reg_Yield<-(1.2*(fixedZvalue-50)+(fixedZvalue-150)^2/200)
    }
    
    if(Yield_Error=="Normal"){
      if(VarFunc=="Const"){
        sig<-10
        
      }
      if(VarFunc=="NonConst"){
        sig<-20/2*fixedZvalue^.2
        
      }  
      Estimates[i,,1]<-dnorm(thisW,reg_Yield,sig)
    }
    
    if(Yield_Error=="Beta"){
      if(VarFunc=="Const"){
        sig<-50
        
      }
      if(VarFunc=="NonConst"){
        sig<-(50*myAPH^.2)
        
      }  
      Estimates[i,,1]<-dbeta(alphaVal/(alphaVal+betaVal)+reg_Yield/sig,alphaVal,betaVal)/sig
      
    }
    
  }

  }

 return(list(Estimates=Estimates,Variance=VarEst)) 
}





