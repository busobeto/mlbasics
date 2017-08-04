# graphSimple.R
#Data Visualizations with R
# 
# Author: PatriciaHoffman
###############################################################################
#http://www.slideshare.net/dataspora/a-survey-of-r-graphics#btnNext
#http://rgm2.lab.nig.ac.jp/RGM2/images.php?show=all&pageID=1396
# at bottom of file many parameters are defined

rm(list=ls())
require(graphics)
#  you may want to run the following code once
#library(help="graphics")
#colours() #  List out colors that are available with R
#demo(graphics)  # this is a great demo of basic graphics capabilities
#This will show you several graphics functions
#press the Enter key to mmove to each new graph

example(points)  # example of plotting points options
?par     # defines many plotting parameters
?Devices # help page for devices
#Open device:  Unix/Linux x11(), windows(), mac: quartz()
# files include pdf, png, jpeg,bmp, tiff ...
# Don't forget to close device
#dev.off()    
library(datasets)
library(help="datasets")
## Make plot appear on screen device

# First Simple Example
with(Orange, {
  plot(age, circumference)
})
with(Orange, {
  plot(age, circumference,col = Tree)
  title(main = "Orange Trees", sub = "age of the tree given as days since 1968/12/31")
  mtext("age vs circumference",side = 3, line=0)
  legend("topleft", c("1","2","3","4","5"),col=c(1,2,3,4,5), pch=1,cex=0.7, title="Tree Number", inset = .05)
})

pdf(file = "C:/Dev/workspace/UCSantaCruzSurveyCourse/DataExploration/myplot.pdf") ## Open PDF device; create 'myplot.pdf' in my working directory
## Create plot and send to a file (no plot appears on screen)
with(Orange, {
  plot(age, circumference,col = Tree)
  title(main = "Orange Trees", sub = "age of the tree given as days since 1968/12/31")
  mtext("age vs circumference",side = 3, line=0)
  legend("topleft", c("1","2","3","4","5"),col=c(1,2,3,4,5), pch=1,cex=0.7, title="Tree Number", inset = .05)
})
dev.off() ## Close the PDF file device
## Now you can view the file 'myplot.pdf' on your computer
oldpar <- par(no.readonly=TRUE) #save the values in par
library(chron)
#  the age of the tree can be put into dates with functions in chron
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1), oma = c(2, 0, 2, 0))
par(mfcol = c(2, 3))
with(Orange, {
  plot(dates(age), circumference,col = Tree,axes=FALSE,cex=1.5)
  axis(1,at=dates(age),labels=dates(age))
  axis(2)
  box()
  title(main = "Orange Trees")
  legend("topleft", c("1","2","3","4","5"),col=c(1,2,3,4,5), pch=1, title="Tree Number", inset = .05)#cex=0.4,
  
  plot(age, circumference,col = Tree,xlab="",type="n")
  for(index in 1:5)
  lines(age[Tree==index],circumference[Tree==index], col = Tree[index*6], lty = index)
  points(age,circumference,bg=Tree,col = "black", pch = 21,cex=1.3)
  title(main = "Color Title",
    xlab = "Faint age",
    col.main = "blue", col.lab = gray(.8),
    cex.main = 1.2, cex.lab = 1.0, font.main = 4, font.lab = 3)
  
  pie(circumference[1:7],col = rainbow(7),radius = 1.0,main="pie chart: circumference")  
  boxplot(circumference,dates(age),col=c("blue",'red'),main = "box plot",
     names=c("circumference","age"))
     rug(jitter(circumference), side = 2)
     rug(jitter(dates(age)), side = 4) 
  hist(circumference, main = "histogram",breaks=seq(from=min(circumference),to=max(circumference),length.out = 8))
  barplot(age[1:7],main="bar plot of age[1:7]",names.arg=dates(age[1:7]),col = rainbow(7))
  })


par(oldpar) 
pairs(iris[1:4],main = "iris data pairs plot", pch = 21, bg = c("red", "green", "blue")[unclass(iris$Species)])



############################################

#            Working with Graphs
#                Create synthetic data
#                from Gaussian Distributions

#############################################

?plot                        #list of available plot functions
methods("plot")              #list of available plot methods
?plot.lm                     #


x <- seq(1,20,0.5)
y <- 1 + 2*x + rnorm(39)
z <- x + rnorm(39,sd=0.5)
model <- lm(y~x)
#mfrow = c(1, 2)  # put two plots on same page side by side
#mar  margines for plot 
#oma  outer margins allow extra space for titles at bottom,left,top,right)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), oma = c(2, 0, 2, 0))
myNumber = 7
plot(x,y,main="Scatter Plot - Linear Model")
points(x,z,col="green")
legend("topleft",c("y","z"),col=c("black","green"),pch=1)
abline(model$coef)
#abline(intercept,slope)
plot(x,y,main="Another Plot - Ablines", cex = 2,pch = 23,col="blue",bg="yellow",
     text(5,35,paste("myNumber = ",myNumber)))
abline(model)
abline(h=20,col = 3,lty = "dotted", lwd = 2) #horizontal line at y = 20 color green
abline(v=10,col = 2, lty = 5) #vertical line at x = 10 color red
mtext("Example Plots", cex=2,side = 3, line = 0,outer = TRUE)
mtext("Inner bottom title",side = 1,line = 2)
mtext("Outer Bottom Title", side = 1, outer=TRUE)
# Lots of plots when the following is executed
# the plot for the lm model just constructed
plot(model,which = c(1:6) )  

# illustration of add=TRUE
par(mfrow=c(1,3))
curve(sin,-2*pi,2*pi,col = "red")
curve(cos,-2*pi,2*pi,col = "blue")
curve(sin,-2*pi,2*pi,col = "red")
curve(cos,-2*pi,2*pi,col = "blue",add=TRUE)

# Pause between graphs - requires a click
#   to see the next graph
# 
par(ask=TRUE) #require a click before each graph

## legends with titles at different locations
plot(x, y, type='n')
legend("bottomright", "(x,y)", pch=1, title="bottomright")
legend("bottom", "(x,y)", pch=1, title="bottom")
legend("bottomleft", "(x,y)", pch=1, title="bottomleft")
legend("left", "(x,y)", pch=1, title="left")
legend("topleft", "(x,y)", pch=1, title="topleft, inset = .05",
       inset = .05)
legend("top", "(x,y)", pch=1, title="top")
legend("topright", "(x,y)", pch=1, title="topright, inset = .02",
       inset = .02)
legend("right", "(x,y)", pch=1, title="right")
legend("center", "(x,y)", pch=1, title="center")


###################################################
###            Plot parameters Examined
###################################################

#       Local parameters for individual plots

#pch is character used for plotting;  pch = 19 is a closed circle
#   filled pch 15,16,17,18,open are pch = 0,1,2,3,4
#    pch 21 through 25 have a border and fill option
#    to set the border col = "blue"; to set the fill color bg="yellow"
#lwd is line width
#cex character (or symbol) expansion; cex = 1.2 make dots 20% larger
#ylim and xlim are the limits on the shown axis
#xlab="...";ylab="..."  labels for the x and y axis
#col is color  col = "blue" col.axis, col.lab, col.main, 
#lty is the line type solid, dashed, etc look at lty=1; through lty=6; 
#type is the way lines are plotted
#type = p gives dots, type=l gives the line, type=0 line with dots, type=b line with breaks and dots
# type = c, s, $, h


#       Global parameters for individual plots

# use the par() function
#las: the orientation of the axis labels on the plot
#bg: the background color
#mar: the margin size (bottom, left, top, right)
#oma: the outer margin size (default is 0 for all sides)
#mfrow: number of plots per row, column (plots are filled row-wise)
#mfcol: number of plots per row, column (plots are filled column-wise)

#par(lty=2,pch=19,col="blue") all this will be in effect for all the graphs that follow 
#par() will tell you what you set
#par("lty") tells you what lty is set to

# Plotting Functions
#plot: make a scatterplot, or other type of plot depending on the class of the object being plotted
#lines: add lines to a plot, given a vector x values and a corresponding vector of y values (or a 2-
#   column matrix); this function just connects the dots
#points: add points to a plot
#text: add text labels to a plot using specified x, y coordinates
#title: add annotations to x, y axis labels, title, subtitle, outer margin
#mtext: add arbitrary text to the margins (inner or outer) of the plot
#axis: adding axis ticks/labels

