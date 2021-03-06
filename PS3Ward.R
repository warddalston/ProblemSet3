################################################
### PS3 - Dalston Ward - February 13, 2014   ###
################################################

#Start everything by emptying the workspace.  (don't use this if concurrently working on other projects in R!)
rm(list=ls())
#Second, set the working directory as appropriate.  This makes the code cleaner and helps to keep track of any saved files.
setwd("~/Documents/WashU 2nd Year/Applied Stats Programming/Feb 6/PS3/ProblemSet3/")
#Third, set the seed for replicability
set.seed(1801)
#Fourth, load necessary libraries:
library(plyr)
library(doMC)
library(multicore)
library(foreach)

###### Section A: Sampling Distributions and P-Values #####

### 1. Make a three dimensional array with dim=c(20,5,100) and fill it with random data.

filler <- rpois(20*5*1000,10) #The random data is sampled from a Poisson Distribution
data <- array(filler,dim=c(20,5,1000)) # It is then stuck into an array of appropriate size
rm(filler) #remove the filler data, just to keep the workspace clean

### 2. Make a function to create Y values (for a linear model).  THe Y values shoudl be a linear combination of the X's plus some normally distributed error.  The output should be a 20 by 1000 array. 

Beta <- matrix(c(1,2,0,4,0),ncol=1) # as given by Jacob

#Ycreator: Function for creating Y value from random data and betas

#This function matrix multiplies data by a set of coefficients, and then adds some random noise to the resulting output.  Useful for generating simulated data for modelling exercises. 

#Input: x - a matrix or vector of dim n by k (in the case of a vector, n=1)
#       beta - a vector of length k

#Output: A vector of length n represented simulated outcome values.  

#Author: Dalston G. Ward

Ycreator <- function(x,beta){
  n <- dim(x)[1]
  noise <- rnorm(n)
  Y <- x%*%beta+noise
}

Yvals <- apply(data,3,Ycreator,Beta) #use apply, selecting the 3rd dimension as the margin, to put my function to use on all the different tables in my simulated data 

dim(Yvals) #check to make sure it has the required dimensions (it does!)

### 3. Run 1,000 regressions across all of this simulated data.  Have as the output a 1000 by 6 matrix of estimated regression coefficients. 

#CoefExtract: Function for fitting and extracting coefficients from many models at once

#Function fits a regression but returns only the coefficients.  Useful when one is fitting a large number of regressions of simulated data which is contained in an array, and is interested only in the distribution of the coefficients.  It uses the arguement "i" to index over the elements of the other arguements.  

#Input: i - a sequential vector from 1 to j, where j is the 3rd dimension length of x
#       Y - a matrix of dim n by j of outcome values
#       x - an array of data of dim n by k by j

#Output: a matrix of regression coefficients with dimensions j by k  

#Author: Dalston G. Ward 

CoefExtract <- function(i,Y,x){
  lm(Y[,i]~x[,,i])[[1]] #the first element of the list output of lm is the vector of coefficients! 
}

BetaHats <- laply(1:dim(data)[3],CoefExtract, Yvals,data) #Use laply, with the argument to apply over as the index. (following Jacob's advice to Tommy on Facebook)
dim(BetaHats) #this is a 1000 by 6 matrix, as required.  

### 4. Create a density plot for each of the 6 coefficients (each of which should have been estimated 1,000 times).  What does this distribution represent? 

#CoefDistrPlotter: Function for making multiple Density plots at once

# This function plots the sampling distributions of beta coefficients from OLS models.  It uses the basic density() function to calculate the density.  THe plot size is limited to extend just beyond the minimum and maximum values in the density.  The x-axis is built to have exactly size tick marks, evently spaced between the minimum and maximum values.  The labels are the same for every plot (I couldn't figure out a general way to make the plots say "beta_0", beta_1, etc while still using expression!)

#Input: x - a numeric vector

#Output: A density plot

#Author: Dalston G. Ward 

CoefDistrPlotter <- function(i,x){
  x <- x[,i]
  xmin <- min(density(x)$x) #this and the next line sets up the plot limit and tick marks
  xmax <- max(density(x)$x)
  plot(density(x),
       xlim=c(xmin-abs(.1*xmin),xmax+abs(.1*xmax)), #just above and below the min and max! 
       xaxt="n", #suppress the default axis so that I can make my own
       main=bquote(paste("Sampling distribution of ", beta[.(i-1)])), #I wish I could figure out a way to loop over subscripts witout a bunch of if loops for this.  
       xlab=bquote(hat(beta[.(i-1)])),
       ylab="Kernal density"
       )
  AxisLabels <- round(seq(xmin,xmax,length.out=6),2) #sequence is perfect for axes, as evenly spaces from one number to another 
  axis(1,at=AxisLabels)
  return(NULL)
}

sapply(1:6,CoefDistrPlotter,BetaHats) #apply it down the columns, should print six plots.  

#These distributions represent the sampling distributions of the regression coefficients, beta_0 through beta_5.  

### 5. Alter your code so that you now collect t-statistics for all 1,000 regressions for all siz coefficients.  

#TStatExtract: Function for fitting and extracting T-statistics from many models at once

#Function fits a regression but returns only the coefficients.  Useful when one is fitting a large number of regressions of simulated data which is contained in an array, and is interested only in the distribution of the t-statistics.  It uses the arguement "i" to index over the elements of the other arguements.  

#Input: i - a sequential vector from 1 to j, where j is the 3rd dimension length of x
#       Y - a matrix of dim n by j of outcome values
#       x - an array of data of dim n by k by j

#Output: a matrix of t-statistics with dimensions j by k  

#Author: Dalston G. Ward 

TStatExtract <- function(i,Y,x){
  mod <- lm(Y[,i]~x[,,i]) #fit a model.  We need to get some pieces out.  
  coef <- mod[[1]] #First, extract the coefficients
  SE <- sqrt(diag(vcov(mod))) #Second, the standard errors.  Recall, these are the square roots of the diagonal elements of the variance-covariance matrix.
  Tstats <- coef/SE #Finally, the t-statistics.  Recall, these are simply the betas divided by the standard errors. 
}

Tstats <- laply(1:dim(data)[3],TStatExtract, Yvals,data) #Use laply, with the argument to apply over as the index. (following Jacob's advice to Tommy on Facebook)

dim(Tstats) #this is a 1000 by 6 matrix, as required.  

### 6. For the 1,000 regressions, calculate how many t-statistics are statistically "significant" (p<=.05) for each variable.  (Make sure you use the right degrees of freedom).  Discuss. 

#SignificanceChecker: Function to check how many t-statistics are statistically significant.

#This function calculates the critical values for a two tailed t-test based on user supplied probabilities.  Defaults to p=c(0.025,0.975), or in other words, a test of significance at the .05 level.  It then tells how many test statistics are larger (or smaller) than these critical values in a vector of test statistics. 

#Input: x: a vector of test statistics
#       n: The number of observations in the model used to generate the test statistic
#       k: The number of coefficients in the model used to generate the test statistic
#       p: The significance levels at which to test, given in terms of probablities.  The default is a two tailed test at the .05 level; hence, it is the vector c(0.025,0.975).  The function accepts both single and two tailed tests.

#Output: a scalar value with the number of test statistics which achieve "significance"

#Author: Dalston G. Ward 

SignificanceChecker <- function(x,n,k,p=c(.025,.975)){
  CritValue <- qt(p,n-k-1) #calculate the critical values
  if(length(CritValue)==2){ #for two-tailed tests. 
  SignifTest <- ifelse(CritValue[1] >= x | x >= CritValue[2],TRUE,FALSE) #compare the T-stats to the critical values.  If a statistic is smaller than the lower critical bound or larger than the upper critical bound, it is signifiacnt, and TRUE is recorded.  Otherwise, FALSE is recorded.  
  }
if(length(CritValue)==1 && CritValue < 0){ # for lower tail tests. 
  SignifTest <- ifelse(CritValue >= x,TRUE,FALSE)
  }
if(length(CritValue)==1 && CritValue > 0){ # for upper tail tests. 
  SignifTest <- ifelse(CritValue <= x,TRUE,FALSE)
}
  sum(SignifTest) #just figure out how many trues, and return.  
}

apply(Tstats,2,SignificanceChecker,20,5) #apply it down the columns.  n=20 because there are 20 observations per regression, and k=5 becuase there are 5 coefficients (not counting the intercept!).

# I see that 1000 out of 1000 coefficient estimates are significantly different than zero for variables 1,2, and 4.  In contrast, I see that only 61 and 48 coefficients are significantly different than zero for variables 3 and 5, respectively.  This is intuitive: recall, we know the true data generating process for the outcome variable.  The X values, which are random draws from a Poisson distribution, were multiplied by the Beta vector c(1,2,0,4,0), to generate the outcome variables.  Thus, we should not expect the coefficient on variables 3 and 5 to be significantly different from zero, because we know the true value of the parameter to be zero! In contrast, we should definitely expect variables 1,2, and 4 to be significantly different from zero, becuase the true values of these coefficients are not zero (they are, instead, 1,2, and 4).

# Furthermore, that we see around 5% of coefficients for variables (intercept), 3, and 5 are significantly different from zero is also not surprising.  The definition of a p value is "the probablity of obtaining a test statistic at least as large in absolute value as the observed test statistic, given that the null hypothesis is true" (Gerber and Green, 2012, page 64).  Our critical values lead us to accept as significant only those test statistics which generate a p-value of .05 or less.  However, these p value of .05 or less mean that there is still a 1 in 20 chance of observing such a test statistic when the null hypothesis is true.  In this case, we know that the null hypothesis is true, and thus observing these 5% (or 1 in 20!) test statistics that are significant is merely an artifact of our null-hypothesis testing procedure.  

### 7. Re-run that code in parallel. Using the system.time command, estimate how much time is saved (or not) using the parallel commmand.  

#Per Jacob's instructions on Facebook, we can rerun any section of out code.  I'll re-run my code for extracting the coefficient vector from a regression.  

#For reference, I time my original code:
OriginalTime <- system.time(BetaHats <- laply(1:dim(data)[3],CoefExtract, Yvals,data))


#Now, I run it in parallel
registerDoMC(cores=4) # I have a 4 core machine.  (Finally figured this out using activity monitor!)
ParallelTime <- system.time(BetaHats <- laply(1:dim(data)[3],CoefExtract, Yvals, data, .parallel=TRUE))

#compare the two times
OriginalTime;ParallelTime 

# we see that the user and system times are both higher for the parallel run, but that the elapsed time is actually lower for parallel. MY real time gain estimate is be about .02 seconds.  I re-run these commands several times to get a distribution of times for parallel and non-parallel. 

####### Note: the next two lines take about 3 minutes each to run! ######
OriginalTimeDist <- replicate(200,system.time(laply(1:dim(data)[3],CoefExtract, Yvals,data))[3])

ParallelTimeDist <- replicate(200,system.time(BetaHats <- laply(1:dim(data)[3],CoefExtract, Yvals,data, .parallel=TRUE))[3])

#Compare the average times:
TimeDiffBar <- mean(OriginalTimeDist)-mean(ParallelTimeDist) 
cat("I estimate that parallel processing is about", TimeDiffBar, "seconds faster than normal processing","\n")

#Is this difference significant?  Run a t-test. 
t.test(OriginalTimeDist,ParallelTimeDist) #Both with and without assuming equal variance, it seems that parallel is significantly faster than normal.

#A boxplot helps show these differences.  
boxplot(OriginalTimeDist,ParallelTimeDist,names=c("Normal","Parallel"),main="Box Plot of System Time distributions",ylab="Elapsed Time",horizontal=TRUE) #Not only do we have further evidence that parallel is faster, but it seems that the variance in elapsed times is lower with parallel processing.  

###### Section B: Calculating Fit Statistics ######


### 1. Using the Out of Step dataset, randomly subset the data into two partitions.  Use one partition (your "training set") to build at least three statistical models where incumbent vote share is the dependent variable.  (Your statistical models can be anything you like, but have some fun with it.  R has lots of machine learning packages, etc. )

#a. Begin by reading in the data from Jacob's website. 
StepData <- read.table("incumbents_0.txt",header=TRUE,sep="\t",row.names=1,stringsAsFactors=FALSE) 

#b Casewise delete all observations with a missing value for the outcome variable
StepData <- StepData[!is.na(StepData$voteshare),]

#c. make the partitions as required.

set.seed(1801) #replicability
Training <- sample(1:nrow(StepData),nrow(StepData)/2,replace=FALSE) #pick out half of the observations in the data at random to be in the training set.

TrainingSet <- StepData[Training,] #Makes the training data
TestingSet <- StepData[-Training,] #everything that isn't training is testing. 

nrow(TrainingSet)+nrow(TestingSet) #Check that every observation of the origninal 6562 ended up in one of the two data sets.  Looks good.  
any(identical(rownames(TrainingSet),rownames(TestingSet))) #Further proof of a successful partition: there are no identical row names in the two data sets.  

#c. fit at least three statistical models. 

#First, an OLS with a year fixed effect.  
Mod1 <- lm(voteshare~inparty+population+urban+seniority+chalquality+south,data=TrainingSet)
summary(Mod1)

#Second, a quasibinomial GLM
Mod2 <- glm(voteshare~inparty+population+urban+seniority+chalquality+south,data=TrainingSet,family=quasibinomial)
summary(Mod2)

#Third, an OLS with a logit transformed response

#Logit and inverse logit function

#These two functions allow one to quickly apply and unapply the logit transformation to a vector

#input: a vector
#output: a vector of equal length

#Author: copied directly from Julien J. Faraway in "Linear Models with R", Chapman & Hall/CRC, 2005. 
logit <- function(x) log(x/(1-x))
ilogit <- function(x) exp(x)/(1+exp(x))

#Now that the funcitons are appropriately defined and documented, run the model.  
Mod3 <- lm(logit(voteshare)~inparty+population+urban+seniority+chalquality+south,data=TrainingSet)
summary(Mod3)

#d. make predictions and the naive forecast

#First, predictions from my OLS
PredMod1 <- predict(Mod1,newdata=TestingSet)

#Second, predictions from my quasibinomial GLM
PredMod2 <- predict(Mod2,newdata=TestingSet,type="response")

#Thid, predictions from my logit transformed OLS
PredMod3 <- ilogit(predict(Mod3,newdata=TestingSet))

#Fourth, make the vector of naive forecasts 
ModNaive <- lm(voteshare~1,data=TrainingSet) #the naive model, based on the training data
PredNaive <- predict(ModNaive,newdata=TestingSet)

#make sure the predictions worked as expected. 
lapply(list(PredMod1,PredMod2,PredMod3,PredNaive),length) #all three prediction vectors are the same length, which is good. 
lapply(list(PredMod1,PredMod2,PredMod3,PredNaive),function(x) sum(is.na(x))) # There are the same number of NA's in each vector of predictions, except for the Naive predictions, which have no missing values.  Will account for this below.    
lapply(list(PredMod1,PredMod2,PredMod3,PredNaive),function(x) which(is.na(x))) # The NA's are the same observations in the three models, which is good.

### 2. Write a funciton that takes as arguments (1) a vector of "true observed outcomes (y), (2) a matrix of predictions (P), and a vector of naive forecasts (r).  The Matrix should be organized so that each column represents a single forecasting model and the rows correspond with each observation being predicted. The function should output a matrix where each column corresponds with one of the above fit statistics and each row corresponds to a model.  

# I begin this problem by defining error calculting functions and fit statistic calculating functions 


# AbsError: Function to calculate the absolute error of a prediciton

#This function calculates the absolute errors of a model's predictions in comparison to observed y values.  It does this according to the formula e = |p-y| where e is the error, p is the prediction, and y is observed.  

#input: p a vector of predictions of length n
#       y a vector of observed values of length n
#       na.rm=TRUE: a logicial, indicating whether or not NA's in the prediction vector should be included in the final output or should be dropped. 

#output: a vector of errors of length n minus any NA's in p and y when N (or length n when na.rm=FALSE)

#Author: Dalston G. Ward 

AbsError <- function(p,y,na.rm=TRUE){
  if(na.rm==TRUE){
    p2 <- p[!is.na(p)] #subset out the missing values from both the predictions and observed values
    y2 <- y[!is.na(p)]
    abs(p2-y2) # and calculate the error. 
  } else {
    #This little function calculates the errors or reports an NA as appropriate. 
    ErrorWithNaCalculator <- function(x,y){ 
      ifelse(is.na(x) | is.na(y),NA,abs(x-y))
    }
  mapply(ErrorWithNaCalculator,p,y) #Applies the above function to the data as appropriate.  
  }
}

#AbsPercentError: Function to calculate the absolute percentage error

#This function calculates the absolute percentage error for a model's predictions in relation to the observed y values.  It does this according to the formula a= |p-y|/|y|*100, where a is the absolute percentage error, p the predictions, and y the observed values.

#input: p a vector of predictions of length n
#       y a vector of observed values of length n
#       na.rm=TRUE: a logicial, indicating whether or not NA's in the prediction vector should be included in the final output or should be dropped. 

#output: a vector of errors of length n minus any NA's in p and y when N (or length n when na.rm=FALSE)

#Author: Dalston G. Ward 

AbsPercentError <- function(p,y,na.rm=TRUE){
  if(na.rm==TRUE){
    e <- AbsError(p,y,na.rm=TRUE) #make use of the function written above to calculate absolute errors, e (and to subset out all of those missing values) 
    e/abs(y[!is.na(p)])*100 
  } else {
    e <- AbsError(p,y,na.rm=FALSE)
    
    #This little function calculates the errors or reports an NA as appropriate. 
    ErrorWithNaCalculator <- function(x,y){ 
      ifelse(is.na(x) | is.na(y),NA,e/abs(y)*100)
    }
    mapply(ErrorWithNaCalculator,e,y) #Applies the above function to the data as appropriate.  
  }
}

#RMSECalc: Function to calculate the RMSE (root mean squared error)

#This function calculates the root mean squared error of a set of predictions in comparison to the true values of y.  It does this according to the formuala sqrt(sum(e)^2/n) where e is the absolute error and n is the number of observations.  

#input: e a vector of absolute errors

#output: a scalar value

#Author: Dalston G. Ward

RMSECalc <- function(e){
  sqrt(sum(e^2)/length(e))
}

#RMSLECalc: Function to calculate the RMSLE (root mean squared log error)

##This function calculates the root mean squared error of a set of predictions in comparison to the true values of y.  It does this according to the formuala sqrt(sum((log(p+1)-log(y+1))^2)/n) where p is the predicted value, y is the observed outcome, and n is the number of observations.  

#input: p a vector of predicted values 
#       y a vector of observed outcomes 

#output: a scalar value

#Author: Dalston G. Ward

RMSLECalc <- function(p,y){
  p2 <- p[!is.na(p)] #subset out the missing values from both the predictions and observed values
  y2 <- y[!is.na(p)]
  LogError <- abs(log(p2+1)-log(y2+1)) # and calculate the error.
  sqrt(sum(LogError^2)/length(LogError))
}

#MRAECalc: Function to calculate MRAE (median relative absolute error)

#This function calculates the median relative absolute error for a set of predictions from a model relative to the errors from the naive predictions.  It does this by calculating median of the relative error, given by the formula e/b, where e is the absolute error, and b is the naive error.

#input: e: a vector of absolute errors
#       b: a vector of naive errors

#output: a scalar value

#Author: Dalston G. Ward

MRAECalc <- function(e,b,P){
  MissingPrediction <- which(!complete.cases(P))
  REs <- e/b[-MissingPrediction]
  median(REs)
}

#FitStatistics: This function calculates several fit statistics for multiple models at once.  

#This function allows for the calculation of several measures of fit for several sets of predictions.  The user must give as inputs a vector of observed y values and a matrix of predictions.  The number of rows in the prediction matrix and the length of the y values should be the same.  It is acceptable to have NA's in the vector of predictions, however.  It can also optionally take a vector of naive predictions, which are used in the calculation of one fit statistic.  The function calculates 6 different fit statistics: RMSE, MAD, RMSLE, MAPE, MEAPE, and MRAE.  The user can specify which of the 6 statistics to calculate; the defualt option is all six.  Additionally, the user can opt to not supply the vector of naive predictions, in which case the MRAE will not be calculated, along with the fit statistics for the naive preditions. 

#Input: y - a vector of observed outcomes
#       P - an n by k matrix of predictions, where n is the length of y and k is the number of models used to generate predictions
#       r - an optional vector of length n with naive predictions.  Defaults to NULL
#       statistic - a character vector specifying which fit statistics the function should calculate.  Defaults to all fit statistics.  The options are RMSE, MAD, RMSLE, MAPE, MEAPE, and MRAE.  MRAE can only be calculated when r is not set to null

#Output: the output is a matrix of fit statistics.  The number of rows is the number of models for which fit statistics are calculated, and the number of columns is the number of fit statistics calculated.  The row names correspond to model names and the column names to fit statistics.  

#Author: Dalston G. Ward 

FitStatistics <- function(y,P,r=NULL,statistic=c("RMSE","MAD","RMSLE","MAPE","MEAPE","MRAE")){
  #The first few lines of the function create objects with the absolute errors, absolute percentage errors, and baseline errors (when these are necesssary)
  E <- apply(P,2,AbsError,y,na.rm=TRUE)
  A <- apply(P,2,AbsPercentError,y,na.rm=TRUE) 
  
  #Create an empty object for each fit statistic.  This will be useful later in putting together the output
  RMSE <- MAD <- RMSLE <- MAPE <- MEAPE <- MRAE <- NULL
    
  #Start calculating the fit statistics
  if(is.null(r)){ #when r, the baseline predictions, are not given it calculates only the first five statistics (further subsetted by what is given as input in the "statistic" arguement).  
  if("RMSE"%in%statistic){ #little if loops allow for the user to specify exactly which statistics to calculate.
    RMSE <- apply(E,2,RMSECalc) #all of the statistics are calculted for multiple sets of predictions at once using apply.  
  }
  if("MAD"%in%statistic){
  MAD <- apply(E,2,median)
  }
  if("RMSLE"%in%statistic){
  RMSLE <- apply(P,2,RMSLECalc,y)
  }
  if("MAPE"%in%statistic){
  MAPE <- apply(A,2,mean)
  }
  if("MEAPE"%in%statistic){
  MEAPE <- apply(A,2,median)
  }
  } else { #From here down only comes into play when r is not null!
    MissingPrediction <- which(!complete.cases(P)) #this helps make the vector of baseline errors the same length as the more sophisticated predictions
    b <- AbsError(r,y) #make the errors for baseline predictions. 
    #The next three lines combine the two types of errors into single objects
    E <- cbind(E,b[-MissingPrediction]) 
    ba <- AbsPercentError(r,y,na.rm=TRUE) 
    A <- cbind(A,ba[-MissingPrediction])
    
    if("RMSE"%in%statistic){ #if loops and apply are used to calculate statistics, as above
      RMSE <- apply(E,2,RMSECalc)
    }
    if("MAD"%in%statistic){
      MAD <- apply(E,2,median)
    }
    if("RMSLE"%in%statistic){
      RMSLE <- c(apply(P,2,RMSLECalc,y),RMSLECalc(b,y)) #the RMSLECalc's input is slightly different from the others, so it must be calculated slightly differently. 
    }
    if("MAPE"%in%statistic){
      MAPE <- apply(A,2,mean)
    }
    if("MEAPE"%in%statistic){
      MEAPE <- apply(A,2,median)
    }
    if("MRAE"%in%statistic){
      MRAE <- apply(E,2,MRAECalc,b,P)
    }
  }

  #this is the output. If a statistic isn't calculated, then it is still null, and as such, doesn't appear in this output matrix
  output <- cbind(RMSE,MAD,RMSLE,MAPE,MEAPE,MRAE)
  if(is.null(r)){ #get the row names right based on whether or not the naive model is considered
  row.names(output) <- colnames(P)
  } else {
  row.names(output) <- c(colnames(P),"ModNaive")
  }
  return(output)
}

#Create the necessary inputs:
Y <- TestingSet[,"voteshare"]
P <- cbind(PredMod1,PredMod2,PredMod3)
r <- PredNaive

FitStatistics(Y,P,r) #give it a whirl.  

###3. Ensure theuser can choose which statistics to calculate and that it works without the null model

### These features are already included.  Here are some examples:

#no naive predictions
FitStatistics(Y,P)

#Only RMSE and MAD
FitStatistics(Y,P,statistic=c("RMSE","MAD"))

#Only MAPE, but with a baseline model
FitStatistics(Y,P,r,statistic="MAPE")

###4. Evaluate the accuracy of the models you fit above using the test set.  

#  Recall, I fit three models using the training data above: an ordinary OLS (Mod1), a quasibinomial GLM (Mod2), and a logit transformed OLS (Mod3).  In problems 2 and 3, above, I used the testing data to make predictions from my models, and then calculated some fit statistics for these predictions.  Now, I turn to the question: how well do my models predict the out of sample data? 

#  First, consider the RMSE statistic.  This measure is in effect the standard deviation of the predictions from the observed values.  We see that all three models produce similar RMSE, around 0.09.  However, the lowest RMSE belongs to the basic OLS. 

#  Second, consider the MAD.  This measures the median absolute deviation.  The statistic doesn't account for predictions that are extremely accurate or extremely inaccurate, something that can be both good and bad.  We see that the GKM model has the smallest MAD, followed by the transformed OLS, and the regular OLS.  

# Third, consider the RMSLE.  This is just like the RMSE expect that errors are logged.  One might want this statistic when non-linearities in prediction errors are suspected.  We see the same pattern as with RMSE: the OLS model has the lowest RMSLE, followed by the GLM and the transformed OLS.  

#  Fourth, consider the MAPE.  This statistic captures the average absolute percentage error, meaning an error of .01 on a value of 1 and an error of .01 on a value of 10 are no longer treated equally.  Again, we see that the OLS model does best, followed by the GLM and the transformed OLS.  

# Fifth, consider the MEAPE.  This is the percentage error equivalent to MAD.  As with MAD, it disregards much of the information in the distribution of errors, for better or worse.  We see that the GLM performs best, followed by the regular OLS and the transformed OLS.  

# Finally, consider the MRAE.  This is the Median relative absolute error.  This measures how much predictive accuracy the model gives us relative to the naive prediction.  As with other median based statistics, it is less sensitive to extreme values, for better or worse. A lower statistic is better here: it means that the amount of predictive power the model gives compared to naive predictions is higher. In contrast with the other statistics, we see that the transformed OLS does best, followed by the normal OLS and the GLM. 

#  We also see that the Naive model performs worse than the other three models on all 6 statistics.  This is reassuring, and tells us that the models fitted above have given us some predictive power. 

#  The fit statistics indicate, in aggregate, that the OLS model seems the best suited to making predictions for data in testing set.  The GLM model seems to be second best, followed by the transformed OLS model, which while still an improvement over naive predictions, is the worst model for predictions.  That said, the fit statistics for all three models are very similar, indicating that there is likely no substantially interesting difference in the predictions made by these models.  However, given that the basic OLS has performed the best, a possible conclusion to draw is that there may be occasions in which the simplest modeling approach may actually be the best.  However, this will certainly depend on the purpose of the model being fitted, the data, and the theory behind the model.  


