# IMPORTANT
# set working directory to the folder containing Bday_problem.R:
source('Bday_problem.R')

# probability of birth on any day of year is uniform
# min number of people (23) to put in a room in order to make
# it more probable for at least one share of Bday 
summary.bday.simulation(no.same.day.bday, number.people = 23, iterations = 1e4)

# non-uniform probs of birth: one day of the week significantly more probable
summary.bday.simulation(no.same.day.bday,birthday.probs = c(2/3,rep(1/18,6)), number.people = 5, iterations = 1e4, boot.iter = 1e3)

# triple birthday now, uniform probs
no.triple <- function(x){
  return(no.same.day.bday(x,3))
}
summary.bday.simulation(no.triple, number.people = 50, iterations = 1e4)

# and quadruple
no.quadruple <- function(x){
  return(no.same.day.bday(x,4))
}
summary.bday.simulation(no.quadruple, number.people = 187, iterations = 1e4)

# EDIT I actually realised that case (50,3,365) does not give ~50%
# hance of occurence and a number which does will probably give us 
# bound which would be quite useless as well

# 3 cases giving ~50% of non-occurence for k=2,3,4. First two
# have similar n^(k+1)/d^k ratios. The last one has a hundred times
# inflated one:
tv.bound.no.occurence(23,2,365)

tv.bound.no.occurence(50,3,365)

tv.bound.no.occurence(187,4,365)

# ratio 1
23^3/365^2
# ratio 2
50^4/365^3
# ratio 3
187^5/365^4

# this might be the reason why the bound for the third case is 
# so hopeless

# let us fix ratio to n^(k+1)/d^k to something in the vicinity of 0.1
# and let's increase n and d simultaneously
k <- 2
d <- floor(exp(2:20))
n <- (0.1*d^k)^(1/(1+k))
n[19]
# then the bounds decrease quite nicely
mapply(function(n,d) tv.bound.no.occurence(n,k,d),n,d)




