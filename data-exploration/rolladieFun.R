#    First Example of User Defined Function
#      Example:  rolladie 
# 
# Author: PatriciaHoffman
###############################################################################
#Parameters appear in procedure definitions; 
#   arguments appear in procedure calls.

# the function rolladie accepts two arguments
#  num.sides (default = 6) indicates 
#    the number of sides of the die
#  num.rolls (default = 1) indicates
#    the number of times the dice is rolled
# rolladie returns the numbers rolled

#Define the Function
rolladie = function (num.sides =6, num.rolls = 1)
{
	simulated.rolls <- trunc(runif(num.rolls)*num.sides+1)
	return(simulated.rolls)
}

#Example Calls
#Call the function
#Notice all the different ways the function can be called
rolladie()
rolladie(num.sides=12)
rolladie(num.rolls=10)
rolladie(num.sides = 12, num.rolls = 10)
rolladie(2,10)

answer <- rolladie()
answer
class(answer)
answer <- rolladie(num.sides = 12, num.rolls = 10)
answer
class(answer)
