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

###### Section A: Sampling Distributions and P-Values #####

### 1. Make a three dimensional array with dim=c(20,5,100) and fill it with random data.

filler <- rpois(20*5*1000,10) #The random data is sampled from a Poisson Distribution
data <- array(filler,dim=c(20,5,1000)) # It is then stuck into an array of appropriate size
rm(filler) #remove the filler data, just to keep the workspace clean

### 2. Make a function to create Y values (for a linear model).  THe Y values shoudl be a linear combination of the X's plus some normally distributed error.  The output should be a 20 by 1000 array. 

Beta <- matrix(c(1,2,0,4,0),ncol=1) # as given by Jacob

# Function for creating Y value from random data and betas

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

#Function for fitting and extracting coefficients from many models at once

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
colnames(BetaHats) <- as.character(0:(ncol(BetaHats)-1))
dim(BetaHats) #this is a 1000 by 6 matrix, as required.  

### 4. Create a density plot for each of the 6 coefficients (each of which should have been estimated 1,000 times).  What does this distribution represent? 

#Function for making multiple Density plots at once

# This function plots the sampling distributions of beta coefficients from OLS models.  It uses the basic density() function to calculate the density.  THe plot size is limited to extend just beyond the minimum and maximum values in the density.  The x-axis is built to have exactly size tick marks, evently spaced between the minimum and maximum values.  The labels are the same for every plot (I couldn't figure out a general way to make the plots say "beta_0", beta_1, etc while still using expression!)

#Input: x - a numeric vector

#Output: A density plot

#Author: Dalston G. Ward 

CoefDistrPlotter <- function(x,i){
  xmin <- min(density(x)$x) #this and the next line sets up the plot limit and tick marks
  xmax <- max(density(x)$x)
  plot(density(x),
       xlim=c(xmin-abs(.1*xmin),xmax+abs(.1*xmax)), #just above and below the min and max! 
       xaxt="n", #suppress the default axis so that I can make my own
       main=expression(paste("Sampling distribution of ", beta)), #I wish I could figure out a way to loop over subscripts witout a bunch of if loops for this.  
       xlab=expression(hat(beta)),
       ylab="Kernal density"
       )
  AxisLabels <- round(seq(xmin,xmax,length.out=6),2) #sequence is perfect for axes, as evenly spaces from one number to another 
  axis(1,at=AxisLabels)
  return(NULL)
}

apply(BetaHats,2,CoefDistrPlotter,) #apply it down the columns, should print size plots.  

#These distributions represent the sampling distributions of the regression coefficients, beta_0 through beta_5.  

### 5. Alter your code so that you now collect t-statistics for all 1,000 regressions for all siz coefficients.  

#Function for fitting and extracting T-statistics from many models at once

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

# Function to check how many t-statistics are statistically significant.

#This function calculates the critical values for a two tailed t-test based on user supplied probabilities.  Defaults to p=c(0.025,0.975), or in other words, a test of significance at the .05 level.  

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

