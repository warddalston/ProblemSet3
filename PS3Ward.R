################################################
### PS3 - Dalston Ward - February 2014       ###
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

this <- lm(Yvals[,1]~data[,,1])
str(this)
class(this[[1]]) #NUMERIC

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

dim(BetaHats) #this is a 1000 by 6 matrix, as required.  


