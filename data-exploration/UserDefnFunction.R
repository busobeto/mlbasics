# Author: PatriciaHoffman
# UserDefnFunction
#
#      At bottom is correlation and covariance
#       CalculateSampleCovariance()
#    Examples include
#      
#      scale() which Uses Prostrate Data
#      apply() uses iris data
#
#      defining an operator - solve Ax = b (A matrix, x & b vector)
#            function bslash() 
#             and operator %!%
#      returning multiple values 
#           multipleReturn()
#      recursive function area()
#           uses one-panel trapezium rule to integrate
#      defining functions for use in plotting 
#            mySillyFunction()
#      variable scope and functions
#             fscope() and fpoint()

#
# Author: PatriciaHoffman
###############################################################################
rm(list=ls())
#myfunction <- function(arg1, arg2, …){
#	statements
#	return(object)
#}
#args can be positional or named
#if return() omitted, results of last statement returned



?source

##################################################

#		Preliminaries: built in trunc function

##################################################
#    ?trunc  (others are ceiling(), floor(), round(x,digits = 0), 
#                                signif(x, digits = 6)
#trunc takes a single numeric argument x and 
#returns a numeric vector containing the integers 
#formed by truncating the values in x toward 0

#     ?runif     runif generates random deviates default - between 0 and 1

#set.seed(pi)

################################
#
#    First Example of User Defined Function
#      Example:  rolladie
#
#################################

#Define the Function
rolladie = function (num.sides =6, num.rolls = 1)
{
	simulated.rolls <- trunc(runif(num.rolls)*num.sides+1)
	return(simulated.rolls)
}

#Call the function
#Notice all the different ways the function can be called
rolladie()
rolladie(num.sides=12)
rolladie(num.rolls=10)
rolladie(num.sides = 12, num.rolls = 10)
rolladie(2,10)



##################################################

#		User Defined Function - scale

##################################################

?scale

# read from a local file
prostate<-read.csv("ProstateCancerDataESL.csv", sep=",",header=T)

# read from a url on the net
myurl <- "http://patriciahoffmanphd.com/resources/data/"
prostate<-read.csv(paste(myurl,"ProstateCancerDataESL.csv",sep=""), sep=",",header=T)

prostate[1:3,]

names(prostate)
dim(prostate)

isTRUE(all.equal(prostate$age,prostate[,3]))

#put predictors into X
X <- prostate[,1:8]
class(X)


prostateScale <- scale(X, center = TRUE, scale = TRUE)
prostateScale[1:3,]

userScale <- function(x, numColm) {
#demean x and bring to unit variance
	for(i in 1:numColm){
		m <- sum(x[,i])
		m <- m/length(x[,i])
		x[,i] <- x[,i]-m
		v <- var(x[,i])
		if(v < 0.0000001) v <- 0.0000001    # include for case v = 0
		x[,i] <- x[,i]/sqrt(v)	# don't want a zero divide
	}
	return(x)
}

prostateUser <- userScale(X,8)
prostateUser[1:3,]

prostateUser <- as.data.frame(prostateUser)
prostateScale <- as.data.frame(prostateScale)


isTRUE(all.equal(prostateUser,prostateScale))

##################################################

#		Example of the apply() function  
#            and ...

##################################################

# using apply()
#apply(X, MARGIN, FUN, ...)  (the ... are the arguments to FUN)
#Returns a vector or array or list of values obtained 
#by applying a function to margins of an array or matrix. 
#
#mean(x, trim = 0, na.rm = FALSE, ...)
#trim the fraction (0 to 0.5) of observations to be trimmed 
#from each end of x before the mean is computed. 
#Values of trim outside that range are taken as the nearest endpoint. 

apply(iris[,1:4], 1, mean) # row means
apply(iris[,1:4], 2, mean) # column means

mydata <- matrix(rnorm(30), nrow=6)
mydata
apply(mydata, 1, mean) # row means
apply(mydata, 2, mean) # column means
apply(mydata, 2, mean, trim=0.2) # column means with options
myDataSet<-read.csv("myDataSet.txt", sep=",",header=T)
rownames(myDataSet)<- paste(rep("row",11) , seq(1,11))
apply(myDataSet, 1, mean) # row means
apply(myDataSet, 2, mean) # column means
apply(myDataSet[,1:11], 2, mean) # column means without last row
apply(myDataSet, 2, mean, trim=0.2) # column means with options

#
#
##    scale training set 
##    some of data are factors
##    scale test set using training parameters
##    the target value is in the first colm 
##    Don't scale the target
#
#
##testScale is matrix
##1	1	1	3
##0	2	1	3
##0	3	1	3
##1	4	1	3
##1	5	1	3
##scaleData is matrix
##1	1	1	3	T	red	5
##0	2	1	3	F	blue	4
##0	3	1	3	T	green	3
##1	4	1	3	T	purple	2
##1	5	1	3	F	yellow	1
#
#

##################################################

#		User Defined Function - MATLAB backslash function

#      from R Introduction by W. N. Venables and D. M Smith
##################################################


#
#As a second example, consider a function to emulate directly the Matlab backslash command,
#which returns the coefficients of the orthogonal projection of the vector y onto the column
#space of the matrix, X. (This is ordinarily called the least squares estimate of the regression
#coefficients.) This would ordinarily be done with the qr() function; however this is sometimes
#a bit tricky to use directly and it pays to have a simple function such as the following to use it
#safely.
#Thus given a n by 1 vector y and an n by p matrix X then X y is defined as (X^(T) X)^(-1) X^(T) y,
# where (XTX)^(-1) is a generalized inverse of X^(')X.

?qr
?qr.coef


bslash <- function(X, y) {
	X <- qr(X)
	qr.coef(X, y)
}


#put response into Y
y <- prostate[,9]


# After the function bslash is created, it may be used in statements such as
regcoeff <- bslash(X, y)

############################################

#            Example bslash
#            Solve Ax = b

#############################################

a <- c(2,1,5,3)
ainverse <- c(3, -1, -5, 2)

A <- matrix(a, ncol = 2, nrow = 2)
Ainverse <- matrix(ainverse, ncol = 2, nrow = 2)

A
Ainverse
A%*%Ainverse
bslash(A, c(1,0))  #[1]  3 -1
bslash(A, c(0,1))  #[1] -5  2


#################################################
#
#  Create a binary operator
#     %!%
#
#
###############################################

"%!%" <- function(X, y) {
	X <- qr(X)
	qr.coef(X, y)
}

A%!%c(1,0)  #[1]  3 -1
A%!%c(0,1)  #[1] -5  2

?lm

regressionModel <- lm(y~. , X)

#  classical ordinary least squares
#beta <- solve(t(x) %*% x) %*% t(x) %*% y


##################################################

#		Examine a Model Object
#                    

##################################################

## EXAMPLE, ADDRESSING AN OBJECT

summary(regressionModel)

str(regressionModel)

regressionModel$call

regressionModel$model$y


## LITTLE MORE COMPLEX

attr(regressionModel$model, "terms")

attr(attr(regressionModel$model, "terms"), "term.labels")


## MAYBE USEFUL

quantile(regressionModel$residuals)



multipleReturn = function (x =6, y = 1)
{
	variabley <- x+y
	variablex <- 10
	z <- cbind(variablex,variabley)
	return(z)
}
z <- multipleReturn()
z
answerx <- z[1]
answery <- z[2]
answerx;answery


basic.stats <- function(x,more=F) {
	stats <- list()
	
	clean.x <- x[!is.na(x)]
	
	stats$n <- length(x)
	stats$nNAs <- stats$n-length(clean.x)
	
	stats$mean <- mean(clean.x)
	stats$std <- sd(clean.x)
	stats$med <- median(clean.x)
	if (more) {
		stats$skew <- sum(((clean.x-stats$mean)/stats$std)^3)/length(clean.x)
		stats$kurt <- sum(((clean.x-stats$mean)/stats$std)^4)/length(clean.x) - 3
	}
	unlist(stats)
}


basic.stats(c(45,2,4,46,43,65,NA,6,-213,-3,-45))
basic.stats(c(45,2,4,46,43,65,NA,6,-213,-3,-45),more=T)


#basic.stats(c(45,2,4,46,43,65,NA,6,-213,-3,-45))
#n     nNAs     mean      std      med 
#11.00000  1.00000 -5.00000 79.87768  5.00000 
#basic.stats(c(45,2,4,46,43,65,NA,6,-213,-3,-45),more=T)
#n      nNAs      mean       std       med      skew      kurt 
#11.000000  1.000000 -5.000000 79.877684  5.000000 -1.638217  1.708149 



f <- function(x) {
	for(i in 1:10) {
		res <- x*i
		cat(x,'*',i,'=',res,'\n')
	}
}
f(7)


#7 * 1 = 7 
#7 * 2 = 14 
#7 * 3 = 21 
#7 * 4 = 28 
#7 * 5 = 35 
#7 * 6 = 42 
#7 * 7 = 49 
#7 * 8 = 56 
#7 * 9 = 63 
#7 * 10 = 70 

##################################################

#		Example: Recursive Function 
#              Integration using trapizoid rule
#                   
#         This also illustrates that you can 
#                pass a function as a parameter


#      from R Introduction by W. N. Venables and D. M Smith
##################################################


#The integrand is evaluated at the end points of the range and in the middle. If
#the one-panel trapezium rule answer is close enough to the two panel, then the latter is
#returned as the value. Otherwise the same process is recursively applied to each panel. The
#result is an adaptive integration process that concentrates function evaluations in regions
#where the integrand is farthest from linear. There is, however, a heavy overhead, and the
#function is only competitive

area <- function(f, a, b, eps = 1.0e-06, lim = 10) {
	fun1 <- function(f, a, b, fa, fb, a0, eps, lim, fun) {
		## function ‘fun1’ is only visible inside ‘area’
		d <- (a + b)/2
		h <- (b - a)/4
		fd <- f(d)
		a1 <- h * (fa + fd)
		a2 <- h * (fd + fb)
		if(abs(a0 - a1 - a2) < eps || lim == 0)
			return(a1 + a2)
		else {
			return(fun(f, a, d, fa, fd, a1, eps, lim - 1, fun) +
							fun(f, d, b, fd, fb, a2, eps, lim - 1, fun))
		}
	}
	fa <- f(a)
	fb <- f(b)
	a0 <- ((fa + fb) * (b - a))/2
	fun1(f, a, b, fa, fb, a0, eps, lim, fun1)
}
                              
x <- seq(0,pi,.05)
y <- cos(x)
plot(x,y)
abline(0,0)
                             # if called with eps = 1.0e-16, lim = 100
area(cos,0, pi)              # [1] -2.630118e-13
sin(pi)- sin(0)

x <- seq(0,pi,.05)
y <- sin(x)
plot(x,y)

area(sin,0, pi)              #[1] 2
- cos(pi)-(- cos(0))         #[1] 2  


# and finally


sillyFunction <- cos
class(sillyFunction)
methods(class="function")
as.list.function(cos)
head(cos)


sillyFunction <- sum(cos, sin)

mySillyFunction <- function(f,g,x){
	h <- f(x) + abs(3*g(2*x))
	return (h)
}

x <- seq(0,pi,.05)
y <- mySillyFunction(sin,cos,x)
plot(x,y)




############################################
#
#       Example:  scope of variable
#
#      Be careful!
#
#           R first uses the value of the variable
#               found in the environment
#             If the variable is not in the environment
#                then it uses the local value
#           R adheres to a set of rules that are called lexical scope.
#   
#     Details can be found in 
#           Section 3.5   Scope of Variables page 22
#           Section 4.3.4 Function Scope page 28
#           http://cran.r-project.org/doc/manuals/R-lang.pdf
#
###############################################

# x is a formal parameter, w is a local variable and z is a free variable.
#
fscope <- function(x) {
	print(y)
	y <- 2*x
	w <- 3*x
	print(w)
	print(x)
	print(y)
	print(z)
}

x <- 1
y <- 7
z <- 3

rm(w)         # w is NOT in the environment
fscope(x)
#y = [1] 7    #  f picks up y=7 from the environment
#w = [1] 3    #      w is locally defined
#x = [1] 1    #      x is the formal parameter
#y = [1] 2    #      y has been changed locally
#z = [1] 3    #      z is the free and is defined by the environment
 
#y is still 7 in the environment
y   #[1] 7

rm(y)  # remove y from the environment and try it again

# now a call to fscope results in an error
#fscope(x)
#Error in print(y) : object 'y' not found

##Delete x (if it exists)
rm(x)
#Here x is declared within the function’s scope of the function, 
#so it doesn’t exist in the user workspace. 
mean(x=1:10)  # function is defined and returns 5.5
#[1] 5.5

# x is not in the environment 
#x   #Error: object 'x' not found



#Run the same piece of code with using the <- operator:	
# Now x is in the environment
mean(x <- 1:10)
x # [1] 1 2 3 4 5 6 7 8 9 10


# And ... one more point
ax <- 1
fpoint <- function(ax) return(TRUE)
fpoint <- fpoint(ax <- ax + 1)
#[1] TRUE
ax
#[1] 1

#Notice that the value of ax hasn’t changed! 
#In R, the value of ax will only change 
#if the argument is evaluated in the function. 
#This can lead to unpredictable behaviour:

ax <- 1
f <- function(ax) if(runif(1)>0.5) TRUE else ax
ax              #[1] 1
f(ax <- ax+1)
ax              #[1] 2
f(ax <- ax+1)   #[1] TRUE
ax              #[1] 2
f(ax <- ax+1)   #[1] TRUE
ax              #[1] 2
f(ax <- ax+1)   
ax              #[1] 3

#  If f returns TRUE the value of ax stays the same
#  If f doesn't return anything the value of ax is increased by one



# One final Note
# The following function is documented according to 
# Google's style guide

CalculateSampleCovariance <- function(x, y, verbose = TRUE) {
	# Computes the sample covariance between two vectors.
	#
	# Args:
	#   x: One of two vectors whose sample covariance is to be calculated.
	#   y: The other vector. x and y must have the same length, greater than one,
	#      with no missing values.
	#   verbose: If TRUE, prints sample covariance; if not, not. Default is TRUE.
	#
	# Returns:
	#   The sample covariance between x and y.
	n <- length(x)
	# Error handling
	if (n <= 1 || n != length(y)) {
		stop("Arguments x and y have invalid lengths: ",
				length(x), " and ", length(y), ".")
	}
	if (TRUE %in% is.na(x) || TRUE %in% is.na(y)) {
		stop(" Arguments x and y must not have missing values.")
	}
	covariance <- var(x, y)
	if (verbose)
		cat("Covariance = ", round(covariance, 4), ".\n", sep = "")
	return(covariance)
}

oldpar <- par(no.readonly=TRUE) # record current setting
par(ask=TRUE)
head(iris)
pairs(iris[1:4], main = "Iris Pairs Plot", pch = 21, bg = c("red", "green3", "blue")[unclass(iris$Species)])
#par(mfrow=c(3,1))

s1 <- "Covariance of Petal Length with Petal Width"
s2 <- "Setosa: "
s3 <- "Versicolor: "
s4 <- "Virginica: "

CalculateSampleCovariance(iris$Petal.Length[1:50], iris$Petal.Width[1:50])
cs<-cov(iris$Petal.Length[1:50], iris$Petal.Width[1:50])
cs
text1 <- sprintf("%s %s %1.3f",s1,s2,cs)
plot(iris$Petal.Length[1:50], iris$Petal.Width[1:50])
mtext(paste(text1), side = 3, line = 0)


CalculateSampleCovariance(iris$Petal.Length[51:100], iris$Petal.Width[51:100])
cver<-cov(iris$Petal.Length[51:100], iris$Petal.Width[51:100])
cver
text2 <- sprintf("%s %s %1.3f",s1,s3,cver)
plot(iris$Petal.Length[51:100], iris$Petal.Width[51:100])
mtext(paste(text2), side = 3, line = 0)

CalculateSampleCovariance(iris$Petal.Length[101:150], iris$Petal.Width[101:150])
cvir<-cov(iris$Petal.Length[101:150], iris$Petal.Width[101:150])
cvir
text3 <- sprintf("%s %s %1.3f",s1,s4,cvir)
plot(iris$Petal.Length[101:150], iris$Petal.Width[101:150])
mtext(paste(text3), side = 3, line = 0)

CalculateSampleCovariance(iris$Sepal.Length[1:50], iris$Sepal.Width[1:50])
cov(iris$Sepal.Length[1:50], iris$Sepal.Width[1:50])
plot(iris$Sepal.Length[1:50], iris$Sepal.Width[1:50])
par(oldpar) # restore settings


CalculateSampleCovariance(seq(1,10,0.5),5*seq(1,10,0.5))
plot(seq(1,10,0.5),5*seq(1,10,0.5))
mtext("cov = 39.58 ; cor = 1 ;")
cov(seq(1,10,0.5),5*seq(1,10,0.5))
cor(seq(1,10,0.5),5*seq(1,10,0.5))

par(oldpar) # restore settings






