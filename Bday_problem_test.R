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
summary.bday.simulation(no.quadruple, number.people = 186, iterations = 1e4)

