# My First R Lesson
rm(list=ls())
## setwd("C:/Users/PatriciaHoffman/workspaceR/TestDataSets")
## myurl <- "http://patriciahoffmanphd.com/resources/data/" 
## C:\Dev\workspaceR\TestDataSets
## setwd("C:/Dev/workspaceR/TestDataSets")
##setwd("C:/Dev/workspace/UCSantaCruzSurveyCourse/DataExploration")
#Author: PatriciaHoffman
###############################################################################**

#####################################**
#                R Documentation from Cran                      ----
#####################################**

# R Reference Card
#	http://cran.r-project.org/doc/contrib/Short-refcard.pdf
#

# An introduction to R is at the following web site:
#       http://cran.r-project.org/manuals.html
#       http://cran.r-project.org/doc/manuals/R-intro.pdf

# Working with Data in R
#       http://cran.r-project.org/doc/manuals/r-devel/R-data.html

# R Language Definition  (assumes knowledge of R)
#       http://cran.r-project.org/doc/manuals/R-lang.pdf

# Index into packages by topic
#       http://cran.r-project.org/web/views/

# Full Reference Manual
#       http://lib.stat.cmu.edu/R/CRAN/doc/manuals/fullrefman.pdf

# Guide to writing a new R Package
#       http://cran.r-project.org/doc/manuals/R-exts.html

#####################################**
#                Access Help through R functions                ----
#####################################**

#  Here are some useful commands
#  
#  To get help on the read function > ??read or > help.search("read")
#      this returns a list of relevant topics
#  To see examples  > example(is.factor) gives examples of is.factor

#  help(help) gives information on the help command
#  To get help > ?is.factor or > help(is.factor)
#  To get help on an arithmetic operator >       ?'%%'
#  To get help on a logical operator >           help("!")
#  List of operator presidence 
#      http://127.0.0.1:27514/library/base/html/Syntax.html
#     help.search() is another function to try
#  Others to try are apropos(), help.search(),help.start()
#  Also see: http://cran.r-project.org/manuals.html
#  To get online help, check out: help(RSiteSearch) 
#      or   RSiteSearch() or RSiteSearch("canonical correlation")
#http://journal.r-project.org/   is another web site to look at



#More web page help
#http://rseek.org  #from Google
#http://stackoverflow.com/
#http://stats.stackexchange.com/


#Documentation for sparce matrix operations 
#   ?SparseM::Ops.matrix.csr

#####################################**
#                Access External R Files                        ----
#####################################**

# source("commands.R") To execute an external file      
# sink("record.lis")   To divert results to an external file

###############
#                Useful Commands                                -----
###############
# to see the data sets loaded with R
# library(help="datasets")
#more commands to investigate:
#      attach(), with(), history(Inf), savehistory(), loadhistory()
#      save(), load()
#run code in Batch mode
#R CMD BATCH options  infile.R  outfile.Rout


#####################################***
#                Accessing Demos, Vignettes, and help            ---- 
#####################################  ****

#Good functions to use to get help
#help.start()
#help(functionname)
#args(functionname)
#example(functionname)
#help.search("pattern")
#help(package = "packagename")
#data()
#vignette()  # to see a list of vignettes
#vignette("vignettename") # to view a particular vignette
#vignette(package="vignettename") # to view a particular vignette
#####################################**
#                Swirl, Sweave, and Shiny                   ----
##################################### **                      ****
#                interactive tutorial on R                  ----
##################################### **                      ****
#
install.packages("swirl")
library("swirl")
swirl()
#https://github.com/swirldev/swirl_courses
# use excape key to get out
# use info() to see other options 

#              Sweave                                     ----
#install.packages("utils")
#library("utils")
#vignette(package="utils") # Lets you know that the Sweave vignette exists
#vignette("Sweave") # provides the Sweave manual
#Sweave provides a flexible framework for mixing text and R/S code for automatic report generation. The basic idea is to replace the code with its output, such that the final document only contains the text and the output of the statistical analysis: however, the source code can also be included. 
#
#Usage
#Sweave(file, driver = RweaveLatex(),
#		syntax = getOption("SweaveSyntax"), encoding = "", ...)
#
#Stangle(file, driver = Rtangle(),
#		syntax = getOption("SweaveSyntax"), encoding = "", ...)

#              Shiny                                        -----
#   package shiny - web application builder
# http://www.rstudio.com/shiny/
# http://rstudio.github.com/shiny/tutorial/

# tutorials from Revolution Analytics
#http://www.rstudio.com/resources/training/online-learning/

#  utilities demos vinettes                                -----

# to see what is in the utility package
# library(help="utils")

#to check out various demos
#demo() # for attached packages
#
### All available demos:
#demo(package = .packages(all.available = TRUE))
#
### Display a demo, pausing between pages
#
#####################################**

#####################################              ****
#                graphics                           ----
#####################################              *****

#  interact with google graphics through r
# https://code.google.com/p/r-google-analytics/


library(help="graphics")
#This will show you several graphics functions
#press the Enter key to mmove to each new graph
demo(graphics)
# this demo will show you how to create a histogram and a box plot
search()  # produces list of packages loaded
 #####################################                      ****
#                Interesting Web Sites for Learning R        ____
#####################################**

#  Documentation for RStudio
# https://support.rstudio.com/hc/en-us/categories/200035113-Documentation

 #     Free Books on R can be found at
#http://freecomputerbooks.com/langRBooks.html

#     R Style Guides:
#http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html#assignment
#https://github.com/hadley/devtools/wiki/Style

# Contains Links to Manuals, books, Journals, User Groups, Mailing Lists, Conferences, 
# http://www.r-project.org/


#tutorial web page
#http://math.illinoisstate.edu/dhkim/rstuff/rtutor.html

#interesting web pages for learning r
# http://www.RDataMining.com
# http://wiki.stdout.org/rcookbook/
# http://www.statmethods.net/
# http://www.r-bloggers.com/
# http://rwiki.sciviews.org/doku.php
# http://is-r.tumblr.com/archive
# http://www.ibm.com/developerworks/linux/library/l-r1/index.html
# http://www.i-programmer.info/programming/other-languages/1706-a-programmers-guide-to-r.html
# http://www.cookbook-r.com/
# http://en.wikibooks.org/wiki/R_Programming

#index to more R resources
# http://news.idg.no/cw/art.cfm?id=F66B12BB-D13E-94B0-DAA22F5AB01BEFE7

#####################################              ****


#    Directory Navigation       move around like linux      ----


#####################################             ****

mydir <- getwd()
getwd()
if (!file.exists("junk")) {  # does directory junk exist?
  dir.create("junk")
}

setwd("./junk")
write.csv(date(),file="theDate.txt")
dir()
setwd("../")
getwd()
list.files("./junk")

# #  move around like linux
# setwd("C:/Dev/workspaceR")
# getwd() # = [1] "C:/Dev/workspaceR"
# setwd("./TestDataSets")
# getwd() # = [1] "C:/Dev/workspaceR/TestDataSets"
# setwd("C:\\Dev\\workspaceR\\TestDataSets")
# getwd() # = [1] "C:/Dev/workspaceR/TestDataSets"
# 

setwd(mydir)
######################################################**



#to remove all the data from the working directory
#good practice to start with a clean workspace
rm(list=ls())
#args("rm")

#####################################**


#                 Read Data                 ----


#####################################**

#available interfaces into R
#RDMS interfaces (packages: RODBC, RJDBCR, MySQL, ROracle, RPostgreSQL, RSQLite, DBI)
#install.packages("MySQL")


#####################################**

# 
#   Read data in from a file             ----
#

#####################################**
mydir <- getwd()
# first set the directory to where the file is located
#setwd("C:/Users/PatriciaHoffman/workspaceR/TestDataSets")
#setwd("C:/Dev/workspaceR/TestDataSets")
# is it at the directory you expected it to be?
getwd()
setwd(mydir)



# to get help reading 
#  help.search("read") # results in lists of read functions
#  ?gdata::read.xls    # results in doc for reading an excel file


help("read.csv")
# ?read.table and ?url
args("read.csv")

#Read from a url on the web
myurl <- "http://patriciahoffmanphd.com/resources/data/"
datax<-read.csv(paste(myurl,"filter60HzNotchblanks.txt",sep = ""), sep=" ",header=F)
#notice that indexing starts with the number one
datax[1:5];length(datax);class(datax)

#Read from a local file
datax<-read.csv("filter60HzNotchblanks.txt", sep=" ",header=F)
datax[1:5];length(datax);class(datax)
datay<-read.csv("http://patriciahoffmanphd.com/resources/data/filter50HzNotch.txt", sep=",",header=F)

#   data is comma seperated values with NO header
datay<-read.csv("filter50HzNotch.txt", sep=",",header=F)
datay[1:5];length(datay)
mydata <-c(datax,datay)                 # the function c is for column combine
class(datax)
mydata                                  # data is a list with
mydata[1:5];length(mydata);class(mydata)    #  lenght = length(datax) + length(datay)

# the first 351 entries in the list is datax 
#    while the last 351 entries is datay

# Reading from other data sources
# There are ways to read excel files into R
# install.packages("RODBC") #for older versions of excel
# install.packages("xlsx")  #for newer versions of excel
# There are ways to use sql like statements with R
# install.packages("sqldf")

#####################################**

# 
#           Coercion                  -----
#

#####################################**
# Change the dimension of data - Change it into a matrix

#  Use the dim function to change a numerical list into a matrix
#     the matrix is filled by columns first
mydata <- c(1,2,3,4,5,10,20,30,40,50)
mydata
class(mydata)
length(mydata)
dim(mydata)
dim(mydata) <- c(5,2)  # c(number of rows, number of columns)
dim(mydata)
class(mydata)
mydata
#   see what is in the 1st row
mydata[1,]
#   see what is in the 1st column
mydata[,1]
#   see what is in the 1st three rows
mydata[1:3,]
#   see the last four entries of the second column
mydata[2:5,2]


# set the dim attribute
# coerce a vector into a matrix

morejunk <- seq(1:10)
morejunk                           # [1]  1  2  3  4  5  6  7  8  9 10
length(morejunk)                   # [1] 10
typeof(morejunk)                   # [1] "integer"
class(morejunk)                    # [1] "integer"
dim(morejunk)                      # NULL

# by changing the dimension attribute, morejunk becomes a matrix
# the values of morejunk are filled up by columns      
attr(morejunk, "dim")<- c(2,5)         # this fills up by columns
#                       c(number of rows, number of columns)
morejunk
#      [,1] [,2] [,3] [,4] [,5]
#[1,]    1    3    5    7    9
#[2,]    2    4    6    8   10

dim(morejunk)                      #[1] 2 5
length(morejunk)                   #[1] 10
typeof(morejunk)                   #[1] "integer"
class(morejunk)                    #[1] "matrix"

#Write Data to a File
write.csv(mydata,file="dataMyFirstRLesson.txt")


int7 <- as.integer(c(1,2,3,4,5))
is.numeric(int7)
is.integer(int7)
class(int7)
num7 <-as.numeric(int7)
class(num7)
class(int7)
?str()
str(int7)  

int7
num7    #[1] 1 2 3 4 5

num7[8]<- 8
num7    #[1]  1  2  3  4  5 NA NA  8

# Delete all the Missing Values
newnum7 <- na.omit(num7)
newnum7
#[1] 1 2 3 4 5 8
#attr(,"na.action")
#[1] 6 7
#attr(,"class")
#[1] "omit"
length(newnum7)#[1] 6
newnum7[8]     #[1] NA
newnum7[6]     #[1] 8

num7 <- num7[1:5]
num7    #[1] 1 2 3 4 5

# Example of an Ordered Factor

survey <- factor(c("poor","best","better","good","excellent","poor","excellent","best","good" ))
results <-factor(survey, levels = c("poor","good","better","best","excellent"),ordered = TRUE)
survey
results
results <- as.numeric(results)
results

#notice again that the index starts with the number one
results[results>2]


###################################################**
### Indexing                                        -----
###################################################**

x <- c(0,-3,4,-1,45,90,-5)
x > 0
#[1] FALSE FALSE  TRUE FALSE  TRUE  TRUE FALSE

x[x>0]               #[1]  4 45 90

x[x <= -2 | x > 5]   #[1] -3 45 90 -5

x[x > 40 & x < 80]  #[1] 45 

x[c(4,6)]            #[1] -1 90
x[1:3]               #[1]  0 -3  4
y <- c(1,4)          
x[y]                 #[1]  0 -1


x[-1]                #[1] -3  4 -1 45 90 -5
x[-c(4,6)]           #[1]  0 -3  4 45 -5
x[-(1:3)]            #[1] -1 45 90 -5

pH <- c(4.5,5,3.3,8.2,6.3)
names(pH) <- c('birch','pine','honeysuckle','lilac','broccoli')
pH

pH['birch']
pH[c('pine','lilac')]

pH <- c(UltraAcid=0,VeryAcid=4.7,Neutral=7,StronglyAlkaline=8.7,VeryBasic=14)


pH['UltraAcid']
pH[c('VeryAcid','VeryBasic')]

pH
#UltraAcid         VeryAcid          Neutral StronglyAlkaline VeryBasic 
#0.0               4.7               7.0     8.7              14.0

#  dealing with NA       ----
x <- c( 1,   2,  NA,     4,  NA,  5)
y <- c("a", "b", NA,   "d",  NA, "f")
z <- c(NA, "zb", "zc","zd",  NA, NA)
good <- complete.cases(x, y)
good
#[1]  TRUE  TRUE FALSE  TRUE FALSE  TRUE
goodAll <- complete.cases(x,y,z)
goodAll

x[good]
#[1] 1 2 4 5
y[good]
#[1] "a" "b" "d" "f"

newdf <- as.data.frame(cbind(x,y,z))

complete.cases(newdf)
completedf <- newdf[complete.cases(newdf),]
completedf
nrow(completedf)
sum(complete.cases(newdf))

########################################**
#
#                 Investigate Objects         -----
#                     
#############################################**
myurl <- "http://patriciahoffmanphd.com/resources/data/"

matrixletters<-read.csv(paste(myurl,"matrixletters.csv", sep=""), sep=",",header=T)
matrixletters

dim(matrixletters)
attributes(matrixletters)

str(matrixletters)

##  Find text where numerics should be

numericf <- as.numeric(levels(matrixletters[,2]))
which(is.na(numericf))
levels(matrixletters[,2])[which(is.na(numericf))]
mytext <- levels(matrixletters[,2])[which(is.na(numericf))]
# in the second column of matrixletters
#   the text "forty two" appears
rowNum <- which(matrixletters[,2] == mytext)
rowNum
matrixletters[rowNum,2]
###################################################**
### Investigate Output Methods                     ----
###################################################**

#using paste to output information
paste("matrixletters is of class ",class(matrixletters), " with dimension = ",
     dim(matrixletters)[1]," by ",dim(matrixletters)[2] , 
     " and column names: ", names(matrixletters)[1], names(matrixletters)[2], 
	 " and ", names(matrixletters)[3])

#using sprintf function from c
s1 <- "matrixletters is of class "
s2 <- " with dimension = "
n1<-names(matrixletters)[1]
n2<-names(matrixletters)[2] 
n3<-names(matrixletters)[3]
sprintf(" %s %s %s %d by %d and column name: %s %s and %s", s1,class(matrixletters),s2,dim(matrixletters)[1],dim(matrixletters)[2],n1,n2,n3)

is.factor(matrixletters[1:8,1])
is.numeric(matrixletters[1:8,1])

is.factor(matrixletters[1:8,2])  #[1] TRUE
is.numeric(matrixletters[1:8,2]) #[1] FALSE

is.factor(matrixletters[1:8,3])  #[1] TRUE
is.numeric(matrixletters[1:8,3]) #[1] FALSE

#   tell R that - means there is No data for that entry
matrixnumbers  <-read.csv(paste(myurl,"matrixletters.csv",sep = ''), sep=",",header=T, na.strings = "-")
matrixnumbers
#
#is.factor(matrixnumbers[1:8,1])
#is.numeric(matrixnumbers[1:8,1])
#
#is.factor(matrixnumbers[1:8,2])
#is.numeric(matrixnumbers[1:8,2])

is.factor(matrixnumbers[1:8,3])  #[1] FALSE
is.numeric(matrixnumbers[1:8,3]) #[1] TRUE

class(matrixnumbers)
str(matrixnumbers)
levels(matrixnumbers$second.clm)
summary(matrixnumbers)
matrixnumbers[4,2]  

objects()
matrixletters
ls()
rm(matrixletters)
matrixletters
objects()
save.image(file = "MyFirstRLesson.RData")
rm(list=ls()) # this starts fresh ...
# no variables exist now
objects()
matrixletters
load("MyFirstRLesson.RData")
objects()

############################################ ****

#             Working with lists             ----

###############################################****
rm(list=ls())

# vectors can only contain the same type of data

# first look at a numeric vector
a <- c(1.1, 2.2, 3.3, 4.4, 5.5)
class(a)

# if the vector has one text entry all the entries are 
# converted to text entries
myVector <- c(1,2,"text", "moreText")
myVector
class(myVector) 


#Lists are indexed heterogeneous collections of other R objects

#create a list using list() function

myList <- list(1,2,"text", "moreText")
myList
class(myList)

list1 <- list(a, sum(a), cbind(a, a*a, a>3))

#the function str() shows the structure of an R object

length(list1)
str(list1)

class(list1)
list1[1];class(list1[1])

list1[2];class(list1[2])

list1[3];class(list1[3])
#[[1]]
#a        
#[1,] 1.1  1.21 0
#[2,] 2.2  4.84 0
#[3,] 3.3 10.89 1
#[4,] 4.4 19.36 1
#[5,] 5.5 30.25 1

#access elements of list

list1[[1]]        #[1] 1.1 2.2 3.3 4.4 5.5
list1[[1]][2]     #[1] 2.2

list1[[2]]        #[1] 16.5

list1[[3]]
list1[[3]][3,2]   #10.89 

#deleting an element of a list

list1[[1]] <- NULL
str(list1)

list1[[1]]

#adding elements to a list
list1 <- list(NULL)
str(list1)

#add elements 1 at a time
list1[[1]] <- a
str(list1)

list1[[2]] <- sum(a)
str(list1)

###################################################**
### Even More Lists                               -----
###################################################**

g <- "My First List"
h <- c(25, 26, 18, 39)
j <- matrix(1:10, nrow = 5)      # this fills up by columns
k <- c("one", "two", "three")

mylist <- list(title = g, ages = h, j, k)   
mylist

mylist$ages  #[1] 25 26 18 39
mylist$ages[2] #[1] 26

my.lst <- list(stud.id=34453, 
		stud.name="John", 
		stud.marks=c(14.3,12,15,19))

my.lst

my.lst[[1]]   # this picks out the number 34453  [1] 34453
my.lst[[3]]

my.lst[1]     # this starts with $stud.id as the name

class(my.lst[3])   #[1] "list"  
class(my.lst[[3]]) #[1] "numeric"   


mode(my.lst[1])
mode(my.lst[[1]])

my.lst$stud.id   #[1] 34453

names(my.lst)
names(my.lst) <- c('id','name','marks')
my.lst

my.lst$parents.names <- c("Ana","Mike")
my.lst

length(my.lst)

other <- list(age=19,sex='male')
lst <- c(my.lst,other)
lst


unlist(my.lst)
my.lst


############################################**

#             Operations                    ----

#############################################**

rm(list=ls()) # this starts fresh ...
              # no variables exist now

aa<-c(1,2,144)
aa
class(aa)

#Manipulating data - simple operations:


aa+10
length(aa)
aa
bb<-c(2,6,12)
my_data_set<-data.frame(attributeA=aa,attributeB=bb)
my_data_set

#Indexing data - note [row,colm]  

my_data_set[,1]
my_data_set$attributeA
my_data_set[,2]
my_data_set$attributeB

my_data_set[3,2]  # row = 3, colm = 2 
my_data_set$attributeB[3]
my_data_set[1:2,]
my_data_set[c(1,3),]


#Matrix Arithmetic:  (in Matlab this is aa ./ bb)
aa
bb
aa/bb
#
#Summary Statistics
#

sum(aa)
prod(aa)
mean(aa)
x<-aa
sum((x-mean(x))^2)/(length(x)-1)
var(x)
yy<- c(1,2,3,4,5,6)
xx <- aa + yy        # length(aa) = 3 while lenght(yy) is 6 ... to do the addition aa is repeated 
xx                   #[1]   2   4 147   5   7 150
x <- 6*aa
x
mean(my_data_set[,1])
median(my_data_set[,1])
sqrt(var(my_data_set[,1]))
head(my_data_set)
summary(my_data_set)

new_data_set <- cbind(c(1,2,3),c(10,20,30))   # column bind
rbind(c(1,2,3),c(10,20,30))                   # row bind
c(c(1,2,3),c(10,20,30))
new_data_set
var(new_data_set[,1])                         # variance of first column
scaleData <- scale(new_data_set)
scaleData                                     # scales each column ... 
var(scaleData[,1])
#
#Write a file out
#
write.csv(my_data_set,"my_data_set_file.csv")
#
#Get a help file
#
?write.csv

input  <- c(aa,bb)
matrixx1<- matrix(data = input, nrow =2,ncol =3)
matrixx1
matrixxx<- matrix(data = c(aa,bb), nrow =2,ncol =3)
matrixxx
#covariance matrix of aa and bb independent vectors.
matrixx2 = var(matrixxx)    # each column of matrixxx has been replaced by the variance of that column
matrixx2
xxx <- c(6, 2, 7, 3, 8, 4)
yy  #[1] 1  2  3  4  5  6  
# pair wise max and min
pmax(xxx,yy)                    #[1] 6 2 7 4 8 6
pmin(xxx,yy)                    #[1] 1 2 3 3 5 4

############################################**
#             complex number in r           ----
#
#     recall that e^(2pi*i) = 1
#############################################**

sqrt(-1)     # error
sqrt(-1 +0i) #works  [1] 0+1i
eulerE <- (1 + 1/20)^20  # lim n -> inf of (1 + 1/n)^n
eulerE
eulerE <- 2.71828182845904523
eulerE^0.1
eulerE^(2*pi*(0+1i))#[1] 1-0i

#   seq(from= ,to= , by = ,length = )

ccc <- seq(1,1.5,0.1)
ccc    #[1] 1.0 1.1 1.2 1.3 1.4 1.5                    
cc1 <- seq(from = 1,by = 0.1, length = 5)
cc1    #[1] 1.0 1.1 1.2 1.3 1.4

ccc <- seq(1,1.5,0.1)
ccc

cc1 <- seq(from = 1,by = 0.1, length.out = 5)
cc1   #[1] 1.0 1.1 1.2 1.3 1.4

seq( 1,by = 0.1, length = 5)
cc2 <- 1:5
cc2         #[1] 1 2 3 4 5
cc3  <- 1:5 -2
cc3         #[1] -1  0  1  2  3
aa

#   rep 
 rep(aa, times = 3)  #[1]   1   2 144   1   2 144   1   2 144
 
 rep(aa, each = 3)   #[1]   1   1   1   2   2   2 144 144 144
 
 ############################################**
 
#            Matrix Operations               ----
 
 #############################################**
 
 
 # multiplication
 
 Ident <- diag((rep(1,times = 3))) 
 Ident
 Z <- matrix(1:9, ncol = 3, nrow = 3)
 Z
 W <- Z %*% Ident    # this is regular matrix multiplication
 W
 W <- Z * Ident      # for matlab users this is the same as .* or by entry multiplication
 W
 Diag2 <- 2*Ident
 Diag2
 W <- Z %*% Ident
 W

 
W <- Z/2
W
W <- Z/Z
W
W <- Z%/% Z
W


firstList <- c(1,2,3)
secondList <- c(5,6,7)
newList <-c(firstList, secondList)
newList
############################################**

#            Solve Ax = b                   ----

#############################################**

a <- c(2,1,5,3)
ainverse <- c(3, -1, -5, 2)

A <- matrix(a, ncol = 2, nrow = 2)
Ainverse <- matrix(ainverse, ncol = 2, nrow = 2)

A
Ainverse


A %*% Ainverse
Ident2 <- A %*% Ainverse

# the matrix Ainverse is the inverse of the matrix A

bfirst <- c(1,0)
bsecond <- c(0,1)

# 


# The first row of Ainverse = Ainverse %*% t(bfirst)
# The second row of Ainverse = Ainverse %*% t(bsecond)


Ainverse
Ainverse %*% bfirst
Ainverse %*% bsecond

# Given Matrix A and Vector b, solve for vector x in the equation A x = b
#      (inverse A)%*% A x = (inverse A)%*% b   => x = (inverseA) %*% b

?solve           
solve(A,bfirst)
solve(A,bsecond)
solve(A,Ident2)

help.search("transpose")
?base::t

A
t(A)

############################################**

#             Investigate Taking Samples    ----

#############################################**

mydata <- rnorm(2000, mean = 0, sd = 1)  # returns 2000 numbers
class(mydata)
attr(mydata, "dim") <- c(1000,2)  
datadf <- as.data.frame(mydata)
head(datadf)
mean(datadf[,1])      #[1] -0.04393697
sample(datadf,10)     # this give an error as sample needs a vector not a data frame
?sample
# here is how to take a sample of 10 rows from a data frame
my_seq <- seq(1,1000)
class(my_seq)
sam<-sample(my_seq,10,replace=T)  # produces ten samples with replacement from  the data = seq(1,1000) 
?seq
my_seq
sam
my_sample<-datadf[sam,] 
my_sample
# my_sample is ten randomly selected samples from data (with replacement)

repSampMean <- function (mydata,sampSize = 10) {
  #Sample Repeatedly and see what happens
  real_mean_colm1<-mean(mydata[,1])
  real_mean_colm2<-mean(mydata[,2])
  store_diff_colm1<-rep(0,10000)
  store_diff_colm2<-rep(0,10000)
  for (k in 1:10000){
    sam<-sample(seq(1,1000),sampSize,replace=T)
    my_sample<-mydata[sam,]
    store_diff_colm1[k]<-abs(mean(my_sample[,1])-real_mean_colm1)
    store_diff_colm2[k]<-abs(mean(my_sample[,2])-real_mean_colm2)
  }
  s1 <- "Column 1 - the difference between the real mean and the mean of the sample = "
  s2 <- "Column 2 - the difference between the real mean and the mean of the sample = "
  cat(s1,mean(store_diff_colm1),"\n",s2,mean(store_diff_colm2)) 
}

repSampMean(datadf)
#Change Sample size from 10 to 100
#   with the larger sample size the mean is more accurate
repSampMean(data,100)



############################################**

#             Simple Linear Regression Example -----

#############################################**
# There is a pdf which details the formula ie y~x
#cran.r-project.org/web/packages/Formula/vignettes/Formula.pdf


set.seed(pi)             # set the seed for the random number generator; pi ~ 3.141593 
x = seq(1,20,0.5)
y = 1 + 2*x + rnorm(39)  # rnorm(39) returns 39 randomly generated samples 
                         #      from Normal distribution with (mean=0,sd=1)
plot(x,y)
model <- lm(y~x)

abline(0.5779,2.0327)    # draws a line with intercept 0.5779 and slope 2.0327

# plot a line 
#     with y intercept = 0.5779
#     and slope        = 2.0327

model

#Call:
#		lm(formula = y ~ x)
#
#Coefficients:
#		(Intercept)            x  
#0.5779       2.0327  

summary(model)

model$coef
abline(model$coef)

fit <- predict(model, newdata = as.data.frame(x))
fitx9 <- fit[9]
# 10.74147 
guess <- 0.5779 + 2.0327*x[9]
#[1] 10.7414
oneMoreWay <- model$coef%*%(c(1,x[9]))
#10.74147

summary(model)
min(y-fit)
max(y-fit)
median(summary(model)$residuals)   #[1] -0.1406118
median(y-fit)                      #[1] -0.1406118

#  classical ordinary least squares
beta <- solve(t(x) %*% x) %*% t(x) %*% y
beta

summary(model)#Displays detailed results for the fitted model
coefficients(model)#Lists the model parameters (intercept and slopes) for the fitted model
confint(model) #Provides confidence intervals for the model parameters (95% by default)
fitted(model) #Lists the predicted values in a fitted model
residuals(model) #Lists the residual values in a fitted model
anova(model) #Generates an ANOVA table comparing two fitted models
vcov(model) #Lists the covariance matrix for model parameters
AIC(model) #Prints Akaikes Information Criterion
plot(model) #Generates diagnostic plots for evaluating the fit of a model
predict(model) #Uses a fitted modelto predict response values for a new dataset




############################################****

#             Investigate Working with Dates    ----

#############################################****
########################################****
# working with dates - Number of days since 1970-01-01

# Dates are represented by the Date class
# Times are represented by the POSIXct or the POSIXlt class
# Dates are stored internally as the number of days since 1970-01-01
# Tmes are stored internally as the number of seconds since 1970-01-01



# Times are represented using the POSIXct or the POSIXlt class
# 
#   POSIXct is just a very large integer under the hood; it use a useful class when you want to store times in something like a data frame
#   POSIXlt is a list underneath and it stores a bunch of other useful information like the day of the week, day of the year, month, day of the month
# 
# There are a number of generic functions that work on dates and times
# 
#    weekdays: give the day of the week
#    months: give the month name
#    quarters: give the quarter number ("Q1", "Q2", "Q3", or "Q4")
# 



# help.start()
# install.packages("chron")
# help(package="chron")
library(chron)

dts <- dates(c("02/27/92", "02/27/92", "01/14/92",
				"02/28/92", "02/01/92"))
dtsmin <- min(dts)    #[1] 01/14/92

# We can add or subtract scalars (representing days) to dates or
# chron objects:
c(dts[1], dts[1] + 10)
# [1] 02/27/92 03/08/92
dts[1] - 31
# [1] 01/27/92

# We can substract dates which results in a times object that
# represents days between the operands:
dts[1] - dts[3]
# Time in days:
# [1] 44
diffdate <- dts[1] - dts[3]
diffdate
as.numeric(diffdate)


dates <- c("02/27/92", "02/27/92", "01/14/92", "02/28/92", "02/01/92")
numday <- as.Date(dates, "%m/%d/%y")
numday   #[1] "1992-02-27" "1992-02-27" "1992-01-14" "1992-02-28" "1992-02-01"
julian(numday)             # represents the number of days after [1] "1970-01-01"
#julian(m, d, y, origin)
startDate <- as.Date(c("01/01/00"),"%m/%d/%y")
julian(numday,startDate)   # represents the number of days after [1] "2000-01-01"

numday

#Times can be coerced 
# from a character string using the as.POSIXlt or 
# as.POSIXct function.

x <- Sys.time()
x               #[1] "2014-04-14 20:37:54 PDT"

p <- as.POSIXlt(x)
names(unclass(p))
## [1] "sec"   "min"   "hour"  "mday"  "mon"
## [6] "year"  "wday"  "yday"  "isdst"
p$sec          #[1] 54.55104



#You can also use the POSIXct format.

x <- Sys.time()
x  ## Already in 'POSIXct' format
## [1] "2013-01-24 22:04:14 EST"
unclass(x)
## [1] 1359083054
x$sec
## Error: $ operator is invalid for atomic vectors
p <- as.POSIXlt(x)
p$sec
## [1] 14.37

#use strptime function if dates are written in a different format
datestring <- c("April 1, 2014 10:40", "December 25, 2013 7:00")
x <- strptime(datestring, "%B %d, %Y %H:%M")
x
## [1] "2014-04-01 10:40:00" "2013-12-25 07:00:00"
class(x)
## [1] "POSIXct" "POSIXt" 


#####################################**


#       Investigate with, aggregate, order, merge, apply, and tapply  ----


#####################################**
Orange
orange <- as.data.frame(Orange)
orange
head(orange)

# Does not print out the summary
with(orange, {
			summary(age, circumference)
			plot(age, circumference)
			points(orange[1:7,2:3],col = 1)
			points(orange[8:14,2:3],col = 2)
			points(orange[15:21,2:3],col = 3)
			points(orange[22:28,2:3],col = 4)
			points(orange[29:35,2:3],col = 5)
		})

#prints out the summary, but variable mystats is NOT available after with statement
with(orange, {
			mystats <- summary(circumference)
			mystats
		})
mystats  #Error: object 'mystats' not found

#thatstats is available after the with statement has finished  notice the << 
with (orange, {
			thatstats <<- summary(circumference)
			thatstats
		})
thatstats

Orange
orange <- as.data.frame(Orange)
orange
summary(orange)
with(orange, {
			sumOrange <- summary(age, circumference)
			plot(age, circumference)
		})

# age is actually a date 
dates(orange$age) 

#####################################***
#       aggregate and order             ----
#####################################**

attach(mtcars)  # 1974 motor trend car road tests
                #[, 2]  cyl  Number of cylinders 
                #[,10]  gear  Number of forward gears  
myCar <- cbind(mtcars$cyl,mtcars$gear,mtcars$mpg,mtcars$wt)
rownames(myCar)<-rownames(mtcars)
colnames(myCar)<-c("cyl", "gear", "mpg", "wt")
head(myCar)
dim(myCar)
myCarOrdered <- with(as.data.frame(myCar), myCar[order(mpg, wt),])
myCarOrdered

aggdata <-aggregate(myCar, by=list(cyl,gear), FUN=mean, na.rm=TRUE)
aggdata   #  the average mpg of a car with 4 cylinders and 3 gears is 21.5

#####################################**
#       apply function               ----
#####################################**
# web page intro to apply function:
# http://nsaunders.wordpress.com/2010/08/20/a-brief-introduction-to-apply-in-r/
# using apply()

orangemx <- cbind(orange$age,orange$circumference)
apply(orangemx, 1, mean) # row means
apply(myCar, 2, mean)    # column means
apply(myCar, 2, mean, trim=0.2) # column means with options

#####################################**
#       tapply function              ----
#####################################**
#tapply(X, INDEX, FUN = NULL, ..., simplify = TRUE)
#X is the object; INDEX is the list of factors; FUN is the function to be applied
orange
tapply(orange$age,orange$Tree,mean) # the average age is 922.1429 for every tree
tapply(orange$circumference,orange$Tree,mean) # the average circumference varies 
     # the average circumference for the third tree is 94


#####################################**

#                Investigate Data using SQL like statements  ----

#####################################**

# using SQL statements
#install.packages("sqldf")
library(sqldf)
newCar <- as.data.frame(myCar)
answ <- sqldf("select * from newCar  where cyl=6 order by mpg",  row.names=TRUE)
answ


#####################################**

#    Create a covariance correlation chart of the iris data ----          

#####################################**
iris <- as.data.frame(iris)
t1 <- subset(iris, iris$Species == "setosa")
t2 <- subset(iris, iris$Species == "versicolor")
t3 <- subset(iris, iris$Species == "virginica")

cov.df <- data.frame(
		"fac" = factor(1:3), 
		"cov" = c(cov(t1$Petal.Length, t1$Petal.Width),
				cov(t2$Petal.Length, t2$Petal.Width),
				cov(t3$Petal.Length, t3$Petal.Width)))

cor.df <- data.frame(
		"fac" = factor(1:3), 
		"cor" = c(cor(t1$Petal.Length, t1$Petal.Width),
				cor(t2$Petal.Length, t2$Petal.Width),
				cor(t3$Petal.Length, t3$Petal.Width)))

chart <- merge(cov.df, cor.df,by.x = "fac", by.y = "fac")
names(chart) <- c("SPECIES","COVARIANCE","CORRELATION")
row.names(chart)<- c("setosa","versicolor","virginica")
chart


#SPECIES  COVARIANCE CORRELATION
#setosa           1 0.006069388   0.3316300
#versicolor       2 0.073102041   0.7866681
#virginica        3 0.048824490   0.3221082
#cor(iris[1:50,3],iris[1:50,4])
#[1] 0.33163


#one more way
#install.packages("plyr")       ----
library("plyr")
require(plyr)

iris_cov<-ddply(iris,.(Species),summarize,COVARIANCE=cov(Petal.Length, Petal.Width))
iris_cor<-ddply(iris,.(Species),summarize,CORRELATION=cor(Petal.Length, Petal.Width))
merge(iris_cov,iris_cor)

#one more way using tapply
agg.cor <-tapply(1:nrow(iris), iris$Species, FUN=function(x) cor(iris$Petal.Length[x], iris$Petal.Width[x]))
agg.cov <-tapply(1:nrow(iris), iris$Species, FUN=function(x) cov(iris$Petal.Length[x], iris$Petal.Width[x]))
iris.stats <- merge(as.data.frame(as.table(agg.cov)), as.data.frame(as.table(agg.cor)), by=c("Var1"))
colnames(iris.stats) <- c("SPECIES", "COVARIANCE", "CORRELATION")
iris.stats