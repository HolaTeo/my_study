##Mini test of BTPS

##Create x and Z values
x<-c(1:100)
z<-c(1:100)

##get true values
my.dataframe<-expand.grid(x,z)
y<-with(my.dataframe,Var1^2+Var2^2+Var1^3*Var2)+rnorm(nrow(my.dataframe),0,50)
x<-my.dataframe$Var1
z<-my.dataframe$Var2

#get random sample
set.seed(9389)
randomSample<-sample(10000,size=500)


##fit the BTPS spline
test<-fitBTPS(y[randomSample],x[randomSample],z[randomSample],knots=c(4,4))
true1derv<-2*test$newX+3*test$newX^2*test$newZ 
true2derv<-2+6*test$newX*test$newZ


##Show it fits the derivatives in the easy case
plot3d(c(test$newX,x),c(test$newZ,z),c(test$Fit,y),col=c(rep("Black",length(test$newX)),rep("Red",length(y))))
plot3d(c(test$newX,test$newX),c(test$newZ,test$newZ),c(true1derv,test$firstDerv),col=c(rep("Black",length(test$newX)),rep("Red",length(test$newX))))
plot3d(c(test$newX,test$newX),c(test$newZ,test$newZ),c(true2derv,test$secondDerv),col=c(rep("Black",length(test$newX)),rep("Red",length(test$newX))))

myVar<-VarianceEstimator(test,CropData = F)

theSD<-sqrt(myVar$variance)
plot3d(c(test$newX,test$newX,myVar$xValues,myVar$xValues),c(test$newZ,test$newZ,myVar$zValues,myVar$zValues),c(true2derv,test$secondDerv,test$secondDerv+1.96*theSD,test$secondDerv-1.96*theSD),col=c(rep("Black",length(test$newX)),rep("Red",length(test$newX)),rep("Green",length(myVar$xValues)),rep("Green",length(myVar$xValues))))
