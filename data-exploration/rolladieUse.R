# Use the function rolladie
#    in the file rolladieFun
# 
# Author: PatriciaHoffman
###############################################################################
rm(list=ls())
objects()

#obtain access to the rolladie function
myDir <- getwd()
dir <- "C:/Users/PatriciaHoffman/workspace/UCSantaCruzSurveyCourse/DataExploration"
setwd(dir)
source("rolladieFun.R")
setwd(myDir)
# delimiters for directories are dependent on operating systems
#for those using Linux
#  check out the R function list.files()
#  and try the argument pattern = "\\.tped"
objects()


#Call the rolladie function
#Notice all the different ways the function can be called
rolladie()
rolladie(num.sides=12)
rolladie(num.rolls=10)
rolladie(num.sides = 12, num.rolls = 10)
rolladie(2,10)

