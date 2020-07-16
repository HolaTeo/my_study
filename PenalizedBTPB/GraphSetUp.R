####Graphs for BTPB

##inputs: simDataEst- Monte Carlo simulation values of the mean from the simStudyFunc
#simDatVar-Monte Carlo simulation values of the variance from the simStudyFunc
#quantileV-quantile values at which to calculate the MC plots of the mean for the penalized BTPB fit


##Returns two graphs 
##One with the MC plots of the mean for the penalized BTPB fit
##One with the heat map of the coverage rate calculated from monte carlo replicates

myGraphs<-function(simDataEst,simDatVar,quantileV=c(.1,.25,.50,.75,.9),Title=NULL){
  
  #obtain the mean estimates and the true estimates
  myMeans<-apply(simDataEst[,,-1], c(1,2), mean,na.rm=T)
  myEst<-melt(myMeans)
  myTrue<-melt(simDataEst[,,1])

  

  ##Set up the graphic plot (ideal for five graphs)
  par(mfrow=c(2,3),mar = c(2,1.15,.6,-1) + 2)

  ##estimate the 2.5 and 97.5 quantile values from the MC simulation
  roundedZ<-round(myEst[,2],4)
  myZValues<-sort(roundedZ)[round(length(roundedZ)*quantileV)]
  percentile975<-apply(simDataEst[,,-1], c(1,2), quantile,na.rm=T,p=.975)
  percentile25<-apply(simDataEst[,,-1], c(1,2), quantile,na.rm=T,p=.025)
  myEst975<-melt(percentile975)
  myEst25<-melt(percentile25)
  
  
  ##plot the values for each quantile along with the true and mean values
  for(i in 1:length(quantileV)){
    par(mar = c(2,2,5,1))
    maxValue<-max(myEst975[roundedZ==myZValues[i],3])
    minValue<-min(myEst25[roundedZ==myZValues[i],3])
    
    plot(myEst[,1][roundedZ==myZValues[i]]*myZValues[i],myEst[,3][roundedZ==myZValues[i]],main=paste0(quantileV[i]*100,"th Percentile (z=",round(myZValues[i],-1),")"),type='l',xlab="",col="gray50",ylim=c(minValue,maxValue),ylab="Yield",lty=1,lwd=3,cex.axis=1.3,cex.title=1.5, cex.main=1.5, cex.sub=1.5)

    
    lines(myEst25[,1][roundedZ==myZValues[i]]*myZValues[i],myEst25[,3][roundedZ==myZValues[i]],col="gray60",lty=2,lwd=3)
    
    
    lines(myEst975[,1][roundedZ==myZValues[i]]*myZValues[i],myEst975[,3][roundedZ==myZValues[i]],col="gray60",lty=2,lwd=3)
    
    lines(myTrue[,1][roundedZ==myZValues[i]]*myZValues[i],myTrue[,3][roundedZ==myZValues[i]],col="black",lwd=2)
    
  }
  
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  ##add a legend
  legend("center",c("True","MC Mean","2.5 & 97.5\nPercentiles"),col=c("black","gray50","gray60"),lty=c(1,1,2), lwd=c(2,2,2),cex=1.5)
  
  
  
  ##configure the heat map graph
  
  #calculate the MC variance (for comparison)
  myMCVar<-apply(simDataEst[,,-1], c(1,2), var,na.rm=T)

  myInd<-array(,dim = dim(simDataEst[,,-1]))

  ##find which MC intervals include the true value and calculate the "non-coverage"
  for(i in 2:dim(simDataEst)[3]){
    
    myInd[,,i-1]<-simDataEst[,,i]+1.96*sqrt(simDatVar[,,i-1])<simDataEst[,,1]|simDataEst[,,i]-1.96*sqrt(simDatVar[,,i-1])>simDataEst[,,1]
    myIndMC<-simDataEst[,,i]+1.96*sqrt(myMCVar)<simDataEst[,,1]|simDataEst[,,i]-1.96*sqrt(myMCVar)>simDataEst[,,1]
    
    
  }
  
  
  myCoverage<-apply(myInd, c(1,2), mean,na.rm=T)
  
  dimnames(myCoverage)<-dimnames(myMeans)
  
  
  #get into form to use with ggplot
  myCovMelt<-melt(myCoverage)

  var1Quant<-sort(myCovMelt$Var1)[c(1,round(length(myCovMelt$Var1)*seq(0,1,length.out=8),0))]
  var2Quant<-sort(myCovMelt$Var2)[c(1,round(length(myCovMelt$Var2)*seq(0,1,length.out=8),0))]
  
  myCovMeltSimp<-subset(myCovMelt,Var1%in%var1Quant&Var2%in%var2Quant)
  names(myCovMeltSimp)[3]<-"Coverage"
  myCovMeltSimp$Coverage<-(1-myCovMeltSimp$Coverage)*100
  
  
  ##apply to the graph
  b <- c(80,85,90,95,100)
  ggplot(myCovMeltSimp, aes(Var1, Var2)) +
    geom_tile(aes(fill = Coverage)) +
    geom_text(size=7,aes(label = (round(Coverage, 1)))) +
    scale_fill_gradientn(limits=c(80,100),
                         colours=c("navyblue", "blue", "light blue","white", "red"),
                         breaks=b, labels=format(b),name="Coverage\nRate")+theme_bw()+ scale_x_continuous(name="Coverage Level") +
    scale_y_continuous(name="Land Quality (APH)")+ theme(text = element_text(size=16) , legend.title=element_text(size=18) , legend.text=element_text(size=14))#+ ggtitle( paste(Title ,"Coverage Rate"))
  
  
  
  
  
}


##inputs: simDataEst- Monte Carlo simulation values of the mean from the simStudyKernel
#simDatVar-Monte Carlo simulation values of the variance from the simStudyKernel
#quantileV-quantile values at which to calculate the MC plots of the mean for the penalized BTPB fit


##Returns two graphs 
##One with the MC plots of the mean for the penalized BTPB fit
##One with the heat map of the coverage rate calculated from monte carlo replicates


###Kernel Heat Maps
myGraphsKern<-function(simDataEst,simDatVar,quantileV=c(.1,.25,.50,.75,.9),Yield_Mean="Linear",VarFunc="Const",Yield_Error="Normal" ,Title=NULL){


  
  
  ##Create a matrix for the true values
true<-simDataEst[,,1]


##create matrices for the lower and upper CI bounds
LowerMatrix<-simDataEst[,,-1]-1.96*simDatVar
UpperMatrix<-simDataEst[,,-1]+1.96*simDatVar


#file a matrix with the coverage interval indicators
temp<-array(,dim=dim(simDatVar))
for(j in 1:dim(simDatVar)[3]){
  temp[,,j]<-(LowerMatrix[,,j]<true)&(UpperMatrix[,,j]>true)
}

Coverage<-apply(temp,c(1,2),mean)
dimnames(Coverage)<-list(dimnames(simDataEst)[[1]],dimnames(simDataEst)[[2]])

test<-melt(Coverage,id.var="CovRate")



##Standardize the yield function since yield will vary based on w and z and will not give a nice grid

stdDev<-1
if(Yield_Mean=="Linear"){
  myMean<-(-25+1.3*test$Var2)
  
}

if(Yield_Mean=="Quad"){
  myMean<-1.2*(test$Var2-50)+(test$Var2-150)^2/200
    
}


if(VarFunc!="Const"){
  stdDev<-(test$Var2^.2)/2
  
}
endPoints<-c(40,40)
if(Yield_Error=="Beta" ){
  endPoints<-c(37.5,62.5)
  
}

test$StdYield<-round(test$Var1-myMean,0)
test<-test[test$Var1>(myMean-endPoints[1]*stdDev)&test$Var1<(myMean+endPoints[2]*stdDev),]


FinalTest<-ddply(test,.(Var2,round(StdYield,-1)),summarize,covRate=mean(value))
names(FinalTest)[2]<-"StdYield"
#return data so it fits a reasonable range
FinalTest<-subset(FinalTest,Var2%in%round(seq(100,300,by = 20))&StdYield%in%round(seq(-40,40,by=10)))
FinalTest<-subset(FinalTest,abs(round(StdYield,-1))<40)


##create the heat map using ggplot2
b <- c(60,70,80,90,95,100)


myPlot<-ggplot(FinalTest, aes((StdYield), (Var2))) +
  geom_tile(aes(fill = covRate*100))+geom_text(size=6,aes(label = (round(covRate*100, 1))))+  scale_fill_gradientn(limits=c(60,100),
                                                                                                                   colours=c("Gray 50","blueviolet", "blue", "light blue","white", "red"),
                                                                                                                   breaks=b, labels=format(b),name="Coverage\nRate")+theme_bw()+ scale_x_continuous(name="Standardized Yield") +
  scale_y_continuous(name="Land Quality (APH)")+ theme(text = element_text(size=16) , legend.title=element_text(size=18) , legend.text=element_text(size=14))


###create the percentile and mean values for the MC graph
percentile975<-apply(simDataEst[,,-1], c(1,2), quantile,na.rm=T,p=.975)
percentile25<-apply(simDataEst[,,-1], c(1,2), quantile,na.rm=T,p=.025)
myMean<-apply(simDataEst[,,-1], c(1,2), mean,na.rm=T)

MYw<-as.numeric(dimnames(simDataEst)[[1]])
##create the MC graph
par(mfrow=c(2,3),mar = c(2,0,.6,-1) + 2)
par(mar = c(4,4,4,1))
for(i in 1:5){
  value<-c(120,160,200,240,280)[i]
  percentage<-c(.1,.25,.5,.75,.9)[i]
  where_Loc<-c(1:dim(simDataEst)[2])[as.numeric(dimnames(simDataEst)[[2]])==value]
  plot(MYw,percentile975[,where_Loc],col="Gray50",type="l",lwd=3,lty=2,ylim=c(0,max(percentile975)),main=paste0(percentage*100,"th Percentile (z=",round(value,-1),")"),xlab="",ylab="",cex.main=1.8,cex.axis=1.5)
  lines(MYw,myMean[,where_Loc],col="Gray50",lwd=3,lty=1)
  lines(MYw,percentile25[,where_Loc],col="Gray50",lwd=3,lty=2)
  
  lines(MYw,true[,where_Loc],lwd=3,lty=1)
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",c("True","MC Mean","2.5th & 97.5th\nPercentile"),col=c("black","gray50","gray60"),lty=c(1,1,3), lwd=c(2,4,4),cex=1.8)


myPlot  
}







