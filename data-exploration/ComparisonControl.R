# Comparison and Control Examples
# 
# Author: PatriciaHoffman
###############################################################################
rm(list=ls())

#            Control Flow Statements
#for (var in seq) statement
#while (cond) statement
#if (cond) statement
#if (cond) statement1 else statment2
#ifelse(cond, statement1, statement2)
#switch (var, list of statements)  (actually switch is an R function)


?rnorm
###################################################
#
##		Logical Expressions
#
###################################################


x <- c(1,1,2,2)
y <- c(2,2,1,1)
z <-  ifelse((x<y),  TRUE,  FALSE) 
z    #TRUE  TRUE FALSE FALSE

y <- c(2,0,3,1)
w <-  ifelse((x<y),  TRUE,  FALSE)
w    #TRUE FALSE  TRUE FALSE

z;w
!w   #FALSE  TRUE FALSE  TRUE    not 
z&w  #TRUE  FALSE FALSE FALSE    and  - both must be TRUE
z&&w #TRUE                             the first element of the vector is TRUE
z|w  #TRUE  TRUE  TRUE FALSE   either or - at least one is TRUE   
z||w #TRUE                               the first element of the vector is TRUE
xor(z,w)#FALSE  TRUE  TRUE FALSE exclusive or - only one is TRUE

##################################################

#		Comparison Examples

##################################################

# remember TRUE = 1 and FALSE = 0

x <- stats::rnorm(20)
x
x < 1
y <- x < 1
y
x[x > 0]


## INTEGER  

X1 <- 4 - 2

Y1 <- 3 - 1

X1 == Y1



# FLOATING  

X1 <- .4 - .2

Y1 <- .3 - .1

X1 == Y1

?gettextf
?sprintf
sprintf("%1.10f",X1)

sprintf("%1.20f",X1)  #[1] "0.20000000000000001110"


sprintf("%1.10f",Y1)

sprintf("%1.20f",Y1) #[1] "0.19999999999999998335"


x1 <- 0.5 - 0.3
x2 <- 0.3 - 0.1
x1;x2
x1 == x2 # FALSE on most machines
identical(x1,x2) #FALSE for my machine
all.equal(x1,x2) #TRUE  for my machine
identical(all.equal(x1, x2), TRUE) # TRUE everywhere

# inside an if clause use isTRUE(all.equal(x1,x2))

isTRUE(all.equal(x1,x2)) # TRUE everywhere


#z <- c(32:126, 160:255) # range of most 8-bit charsets, Latin-1 in Unicode
#x <- if(l10n_info()$MBCS) {
#			intToUtf8(z, multiple = TRUE)
#		} else rawToChar(as.raw(z), multiple= TRUE)
### by number
#writeLines(strwrap(paste(x, collapse=" "), width = 60))
### by locale collation
#writeLines(strwrap(paste(sort(x), collapse=" "), width = 60))

# individual variables
pi;355/113
all.equal(pi, 355/113)   #[1] "Mean relative difference: 8.491368e-08"
# not precise enough (default tol) > relative error
all.equal(pi, 355/113, tolerance = .Machine$double.eps ^ 0.5) #[1] "Mean relative difference: 8.491368e-08"
all.equal(pi, 355/113, tolerance = 8.4e-08)                   #[1] "Mean relative difference: 8.491368e-08"
all.equal(pi, 355/113, tolerance = 8.5e-08)                   #[1] TRUE

isTRUE(all.equal(pi, 355/113, tolerance = 8.4e-08))           #[1] FALSE
isTRUE(all.equal(pi, 355/113, tolerance = 8.5e-08))           #[1] TRUE

# vectors

# d45 = (1 1/4 pi, 2 1/4 pi, ...10 1/4pi)
# tan(45)= (1,1,1,1,1,1,1,1,1,1)
# rep(1,10 = (1,1,1,1,1,1,1,1,1,1)
d45 <- pi*(1/4 + 1:10)
all (tan(d45) == rep(1,10)) # FALSE, since not exactly  
all.equal(tan(d45), rep(1,10), tol=0) # to see difference  [1] "Mean relative difference: 1.29526e-15"

isTRUE(all.equal(tan(d45),rep(1,10))) # TRUE

(tan(d45) == rep(1,10))               # [1] FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE



a <- 7; b <- 3        # %% is modulo = remainder after interger division  7/3 = 2 remainder 1
isTRUE(all.equal(a, (a %% b) + b * ( a %/% b )))  #TRUE a = remainder + b*int(a/b) 


# comparisons
aa<-(c(3,2,1,NA,NaN))
aa
dim(aa) <- c(5,1)
aa
temp <- aa>2
temp
tempno <- temp[!is.na(temp)]
tempno
tempyes <- temp[is.na(temp)]
tempyes

naIndex <-which(is.na(aa))
naIndex
valueIndex<- which(!is.na(aa))
valueIndex


new <- (rep(0,nrow(aa)))
new
new[naIndex] <- 1
new
bb <- cbind(aa,new)
bb
bbframe <- as.data.frame(bb)
bbframe

?Quotes
labs <- paste(c("X","Y"), 1:10, sep="")
labs

# sprintf and pi
## various formats of pi :

sprintf("%f is the number pi", pi)
sprintf("%.3f", pi)
sprintf("%1.0f", pi)
sprintf("%5.1f", pi)
sprintf("%05.1f", pi)
sprintf("%+f", pi)
sprintf("% f", pi)
sprintf("%-10f", pi) # left justified
sprintf("%e", pi)
sprintf("%E", pi)
sprintf("%g", pi)
sprintf("%g",   1e6 * pi) # -> exponential
sprintf("%.9g", 1e6 * pi) # -> "fixed"
sprintf("%G", 1e-6 * pi)

#
## Example for identical
#
#identical(1, NULL) ## FALSE -- don't try this with ==
#identical(1, 1.) ## TRUE in R (both are stored as doubles)
#identical(1, as.integer(1)) ## FALSE, stored as different types
#x <- 1.0; y <- 0.99999999999
### how to test for object equality allowing for numeric fuzz :
#(E <- all.equal(x,y))
#isTRUE(E) # which is simply defined to just use
#identical(TRUE, E)
### If all.equal thinks the objects are different, it returns a
### character string, and the above expression evaluates to FALSE
### even for unusual R objects :
#identical(.GlobalEnv, environment())
#### ------- Pickyness Flags : -----------------------------
### the infamous example:
#identical(0., -0.) # TRUE, i.e. not differentiated
#identical(0., -0., num.eq = FALSE)
### similar:
#identical(NaN, -NaN) # TRUE
#identical(NaN, -NaN, single.NA=FALSE) # differ on bit-level

##################################################

#		for Loop Example

##################################################

for(i in 1:5) print(1:i)
i
for(n in c(10,100,1000,200000,500000)) {
	y <- stats::rnorm(n)
	cat(n,":",mean(y),var(y),"\n")
}
letters
f = factor(sample(letters[1:5], 10, replace=TRUE))
f
for( i in unique(f) ) print(i)
i
sort(f)
for (myLetter in letters) print(myLetter)

##################################################
#
#		while Loop Example
#           switch (statement, list)
#           
#
##################################################
keepGoing <- TRUE
count <- 0
food <- junk
myFood <- c("fruit", "vegetable", "meat")
while(keepGoing){
	count <- count + 1
	print(count)
	switch(myFood[count],fruit = (food = "apple"), vegetable = (food = "broccoli"), meet = (food = "chicken"), (food = "nothing"))
	print(paste("The", myFood[count], "I ate today is ",  food))
	# put in a hard stop
	if(count>3.0)keepGoing <- FALSE
}
#?switch
#switch(EXPR, ...)
#switch evaluates EXPR and accordingly chooses one of the further arguments (in ...). 
#EXPR an expression evaluating to a number or a character string. 
#... the list of alternatives.  


##################################################

#		if Statement Example

##################################################

mylogic <- TRUE
if(mylogic)print(pi)else print(2*pi)
mylogic <- FALSE
if(mylogic)print(pi)else print(2*pi)


w <- vector(mode = "numeric", length = 21)
x[21] <- 0
x
for(i in 1:21){
	if(x[i]<0) w[i]=0
	if(x[i]>0) w[i]=1	
	if( isTRUE(all.equal(x[i],0, tol = 0.1))) w[i] <- 100
}

x; w; i

##################################################

#		which - example for indexing

##################################################
set.seed(pi)
x <- stats::rnorm(20)
rm(w)
w <- vector(mode = "numeric", length = 21)  # w is a vector with 21 numerical zeros as entries
w

x
which(x>0)
w[which(x>0)]   = 1    # set  w[index] to 1 if  x[index] is greater than zero
w[-which(x>0)]  = 0    # set  w[index] to 0 if  x[index] is not greater than zero
w[which.min(x)] = 100  # set w[index] to 100 for the index in which x[index] is the smallest

w
which.min(x)  # 9
w[9]

mylist <- c(1,2,3,4,5,6,7,NA)
small<-mylist[which(mylist < 4)]
small
large <-mylist[which(mylist >= 4)]
large
# use which to index create
#use isTRUE(all.equal(x1,x2))to compare

##################################################

#		Which Example - divide sample into test and train sets

##################################################


prostate<-read.csv("ProstateCancerDataESL.csv", sep=",",header=T)

dim(prostate)
names(prostate)
attributes(prostate)

# If the last attribute is FALSE then it is in the test set

I <- seq(from = 1, to = nrow(prostate))
Itest <- which(prostate[I,10] == FALSE)
y.train<-prostate[-Itest,9]
junk<-prostate$lpsa[-Itest]

isTRUE(all.equal(y.train,junk))

dim(prostate)

# funky way to add long comment
if(0){
	"add a bunch of text here" 
     }


dim(prostate)