################################################
### PS3 - Dalston Ward - February 2014       ###
################################################

#Start everything by emptying the workspace.  (don't use this if concurrently working on other projects in R!)
rm(list=ls())
#Second, set the working directory as appropriate.  This makes the code cleaner and helps to keep track of any saved files.
setwd("~/Documents/WashU 2nd Year/Applied Stats Programming/Feb 6/PS3/ProblemSet3/")
#Third, set the seed for replicability
set.seed(1801)

# Section A: Sampling Distributions and P-Values

# 1. Make a three dimensional array with dim=c(20,5,100) and fill it with random data.

filler <- rpois(20*5*1000,10) #The random data is sampled from a Poisson Distribution
data <- array(filler,dim=c(20,5,1000)) # It is then stuck into an array of appropriate size
rm(filler) #remove the filler data, just to keep the workspace clean

# Make a function to create Y values (for a linear model).  THe Y values shoudl be a linear combination of the X's plus some normally distributed error.  The output should be a 20 by 1000 array. 

Beta <- matrix(c(1,2,0,4,0),ncol=1) # as given by Jacob

# Y value Generator 