DailyCases<-read.table("Data_Group3.txt",header = TRUE) #read in data
DailyCases$day<-rep(0,98)
counter<-1
for (j in 1:98) #add in days variable
{
  DailyCases$day[j]<-counter
  counter<-counter+1
}
DailyCases$CumulativeCases<-rep(0,98) #add in total number of cases variable
colnames(DailyCases)=c("DailyCases","day","CumulativeCases") #give appropriate names
counter2<-0
for (j in 1:98)
{
  counter2<-counter2+DailyCases$DailyCases[j] #sum of cases up to day j
  DailyCases$CumulativeCases[j]<-counter2
}

#Produce some plots of the daily cases data and their cumulative sum
par(mfrow=c(1,2)) #create 1x2 frame for plots
plot(DailyCases$day,DailyCases$DailyCases,xlab="Day",ylab="Daily cases",main="Daily cases vs. day",cex=0.8,pch=1) #make points slightly smaller
plot(DailyCases$day,DailyCases$CumulativeCases,xlab="Day",ylab="Cumulative cases",main="Total cases vs. day",cex=0.8,pch=1,col="blue")

#Now define the system of ODEs
library(deSolve)
seir.model<-function(t,x,params){
  beta=params[1]; gamma=params[2]; p=params[3] #take these from input vector, quantities of interest
  N<-x[1]+x[2]+x[3]+x[4] #keeps N up-to-date so dead removed from population
  dS<-((-1)*beta*x[1]*x[3])/N #change in susceptibles
  dE<-(beta*x[1]*x[3])/N-0.25*x[2] #change in exposed
  dI<-0.25*x[2]-gamma*x[3] #change in infected
  dR<-(1-p)*gamma*x[3] #change in recovered
  dD<-p*gamma*x[3] #change in dead
  dC<-0.25*x[2] #change in cumulative infected
  list(c(dS,dE,dI,dR,dD,dC)) #return list of all these
}

initial_values=c(999,0,1,0,0,1) #list of initial values (S(0),E(0),I(0),R(0),D(0),C(0))
param.start=list(beta=0.5,gamma=0.25,p=0.005) #initial guesses for parameters (beta,gamma,p)
fills<-rep(0,98) #make a zero vector to fill covidData with, so columns can have same names
covidData<-data.frame(DailyCases$day,fills,fills,fills,fills,fills,DailyCases$CumulativeCases)
colnames(covidData)=c("day","S","E","I","R","D","C") #give appropriate names
library(minpack.lm)
nls.fit=nlsLM(C~rk(initial_values,day,seir.model,c(beta=beta,gamma=gamma,p=p))[,7],
              start=param.start,data=covidData,lower=c(0.0000001,0.0000001,0.0000001),
              upper=c(1,1,1)) #fit a non-linear regression model
#This fit seems to be reasonably good in terms of what we'd expect of the rates
coef(nls.fit) #coefficients of the model
vcov(nls.fit) #variance of these coefficients
#Now fit the model to the data
t<-1:98 #number of days
param.start_fit=c(beta=coef(nls.fit)[1],gamma=coef(nls.fit)[2],p=coef(nls.fit)[3]) #to fit, need vector not list
solution_fit=as.data.frame(rk(initial_values,t,seir.model,param.start_fit))
colnames(solution_fit)=c("day","S","E","I","R","D","C") #give appropriate names
cumulative_fit<-solution_fit$C #get the cumulative cases from the fit
#Now need to 'unzip' the cumulative cases for the daily cases
daily_fit<-rep(0,98) #make a vector to fill with daily cases
for (i in 0:97)
{
  if (i==97)
  {
    daily_fit[1]<-cumulative_fit[1]
  }
  else
  {
    daily_fit[98-i]<-cumulative_fit[98-i]-cumulative_fit[97-i] #find difference
  }
}
par(mfrow=c(1,1)) #reset frame so plot takes up the whole space
plot(DailyCases$day,DailyCases$DailyCases,xlab="Day",ylab="Daily cases",main="Data and fitted line",cex=0.4,pch=1) #make points slightly smaller
lines(t,daily_fit,col="red")
legend('topright',legend=c("Data","Fitted"),col=c("black","red"),pt.lwd=1,lwd=1,bty="n",text.font=.5)
#Okay, so the fit has the right shape but peaks not sharp enough. May be different with different start vectors

#Check here that the errors about normally distributed (can ignore this section, should be considering errors of all parts of the model)
residuals<-DailyCases$DailyCases-daily_fit
hist(residuals,breaks=20) #more negative than desired
plot(t,residuals,main="Residual plot",xlab="Day",ylab="Residual") #not the best
#Histogram for R0
set.seed(31415) #set seed for reproducibility
betaSamples<-rnorm(100000,mean=coef(nls.fit)[1],sd=sqrt(vcov(nls.fit)[1,1])) #sample betas
gammaSamples<-rnorm(100000,mean=coef(nls.fit)[2],sd=sqrt(vcov(nls.fit)[2,2])) #sample gamma
R0Samples<-rep(0,100000) #initialise vector for R0 values
for (i in 1:100000)
{
  R0Samples[i]<-betaSamples[i]/gammaSamples[i]
}
hist(R0Samples,breaks=100,main="Histogram of R0 values",xlab="R0") #plot histogram with 50 breaks
min(R0Samples) #the minimum that is given with 10000 points is 1.321968
length(which(R0Samples<1.3)) #there are 9 values less than 1.5
R0bestEstimateOrder1<-coef(nls.fit)[1]/coef(nls.fit)[2] #the best estimate using first order Taylor expansion is 2.992241

#Confidence intervals by first method using the data
R0Samples<-sort(R0Samples)
Interval<-c(R0Samples[2500],R0Samples[97500])
Interval
#Katz method for confidence interval (nicely lines up with what we have)
logmean<-log(coef(nls.fit)[1]/coef(nls.fit)[2])
logvariance<-vcov(nls.fit)[1,1]/coef(nls.fit)[1]^2+vcov(nls.fit)[2,2]/coef(nls.fit)[2]^2
lowerKatz<-exp(logmean-qnorm(0.975)*sqrt(logvariance))
upperKatz<-exp(logmean+qnorm(0.975)*sqrt(logvariance))
intervalKatz<-c(lowerKatz,upperKatz)
#Geary-Hinkley transformation for confidence interval (even more nicely lines up)
a<-coef(nls.fit)[2]^2-qnorm(0.975)^2*vcov(nls.fit)[2,2]
b<-2*(vcov(nls.fit)[1,2]*sqrt(vcov(nls.fit)[1,1]*vcov(nls.fit)[2,2])*qnorm(0.975)^2-coef(nls.fit)[1]*coef(nls.fit)[2])
c<-coef(nls.fit)[1]^2-qnorm(0.975)^2*vcov(nls.fit)[1,1]
TSolutionUpper<-(-b+sqrt(b^2-4*a*c))/(2*a)
TSolutionLower<-(-b-sqrt(b^2-4*a*c))/(2*a)
a2<-coef(nls.fit)[2]^2+qnorm(0.975)^2*vcov(nls.fit)[2,2]
b2<--2*(vcov(nls.fit)[1,2]*sqrt(vcov(nls.fit)[1,1]*vcov(nls.fit)[2,2])*qnorm(0.975)^2+coef(nls.fit)[1]*coef(nls.fit)[2])
c2<-coef(nls.fit)[1]^2+qnorm(0.975)^2*vcov(nls.fit)[1,1]
determinant<-b2^2-4*a2*c2 #There is no real solution on this side of the interval
c2 #This is positive, can prove by completing the square that all values satisfy inequality
intervalGeary<-c(TSolutionLower,TSolutionUpper)
#All three confidence intervals coincide, see wikipedia/notes for explanations of above


#We get different values of the parameters with different starting vectors
#Let's do a sort of grid search to see what kind of values we can get for coefficients and variance
initialbeta<-seq(0.1,1,by=0.1)
initialgamma<-seq(0.1,1,by=0.1)
searchInitials=matrix(0,nrow=100,ncol=6)
rownumber<-1
for (i in 1:10)
{
  for (j in 1:10)
  {
    initialcoefs<-c(initialbeta[i],initialgamma[j])
    searchInitials[rownumber,1:2]<-initialcoefs
    param.start1=list(beta=initialbeta[i],gamma=initialgamma[j],p=0.005)
    nls.fit1=nlsLM(C~rk(initial_values,day,seir.model,c(beta=beta,gamma=gamma,p=p))[,7],
                   start=param.start1,data=covidData,lower=c(0.0000001,0.0000001,0.0000001),
                   upper=c(1,1,1)) #fit a non-linear regression model
    finalcoefs<-coef(nls.fit1)[1:2]
    searchInitials[rownumber,3:4]<-finalcoefs
    finalvars<-c(vcov(nls.fit1)[1,1],vcov(nls.fit1)[2,2])
    searchInitials[rownumber,5:6]<-finalvars
    rownumber<-rownumber+1
  }
}
#For most of these initial values the fitted value of beta is 1, seems that R0 is going to be about 3
plot(DailyCases$day,DailyCases$DailyCases,xlab="Day",ylab="Daily cases",main="Data/fitted lines",cex=0.4,pch=1,ylim=c(0,50)) #make points slightly smaller
lines(t,daily_fit,col=1)
colour<-2
sensibleRows<-which(searchInitials[,3]<1)
dailyFits=matrix(0,nrow=length(which(searchInitials[,3]<1)),ncol=98)
for (i in 1:length(which(searchInitials[,3]<1)))
{
  row<-sensibleRows[i]
  param.start_fitted=c(beta=searchInitials[row,3],gamma=searchInitials[row,4],p=0.005)
  solution=as.data.frame(rk(initial_values,t,seir.model,param.start_fitted))
  cumulative<-solution$'6'
  daily_fitted<-rep(0,98) #make a vector to fill with daily cases
  for (j in 0:97)
  {
    if (j==97)
    {
      daily_fitted[1]<-cumulative[1]
    }
    else
    {
      daily_fitted[98-j]<-cumulative[98-j]-cumulative[97-j] #find difference
    }
  }
  dailyFits[i,]<-daily_fitted
  lines(t,daily_fitted,col=colour)
  colour<-colour+1
}
#The third row in dailyFits has the biggest peak
sensibleRows[3] #41st column in searchInitials gives the big peak
searchInitials[sensibleRows[3],] #the final fitted parameters are beta=0.739 and gamma=0.164, give larger value for R0
lines(t,dailyFits[3,],col="black") #check that this changes the colour to black (it does)
#Pull estimates and their variances so we can look at R0
meanAndVariances<-searchInitials[sensibleRows[3],3:6]
logmean2<-log(meanAndVariances[1]/meanAndVariances[2])
logvariance2<-meanAndVariances[3]/(meanAndVariances[1]^2)+meanAndVariances[4]/(meanAndVariances[2]^2)
lowerKatz2<-exp(logmean2-qnorm(0.975)*sqrt(logvariance2))
upperKatz2<-exp(logmean2+qnorm(0.975)*sqrt(logvariance2))
intervalKatz2<-c(lowerKatz2,upperKatz2) #Katz interval massively wide
#Note that we can't use the data here because its actually very likely that gamma will be negative
#Need covariances if we want to get Geary interval
param.start2=list(beta=0.5,gamma=0.1,p=0.005)
nls.fit2=nlsLM(C~rk(initial_values,day,seir.model,c(beta=beta,gamma=gamma,p=p))[,7],
               start=param.start2,data=covidData,lower=c(0.0000001,0.0000001,0.0000001),
               upper=c(1,1,1)) #fit a non-linear regression model


aGeary<-coef(nls.fit2)[2]^2-qnorm(0.975)^2*vcov(nls.fit2)[2,2]
bGeary<-2*(vcov(nls.fit2)[1,2]*sqrt(vcov(nls.fit2)[1,1]*vcov(nls.fit2)[2,2])*qnorm(0.975)^2-coef(nls.fit2)[1]*coef(nls.fit2)[2])
cGeary<-coef(nls.fit2)[1]^2-qnorm(0.975)^2*vcov(nls.fit2)[1,1]
TSolutionLowerGeary<-(-bGeary-sqrt(bGeary^2-4*aGeary*cGeary))/(2*aGeary)
#First inequality holds for t>1.37, this is the Geary interval in this case
a2Geary<-coef(nls.fit2)[2]^2+qnorm(0.975)^2*vcov(nls.fit2)[2,2]
b2Geary<--2*(vcov(nls.fit2)[1,2]*sqrt(vcov(nls.fit2)[1,1]*vcov(nls.fit2)[2,2])*qnorm(0.975)^2+coef(nls.fit2)[1]*coef(nls.fit2)[2])
c2Geary<-coef(nls.fit2)[1]^2+qnorm(0.975)^2*vcov(nls.fit2)[1,1]
#Second inequality holds for all t
#This is a bad approximation (see wikipedia for reason why)