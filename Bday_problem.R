require('boot')
seq.alt <- function(a, b){
  # helper function, returns empty sequence if a > b,
  # otherwise just normal sequence
  if (a <= b)
    return(a:b)
  else
    return(NULL)
}

simulate.room<- function(statistics, birthday.probs = rep(1,365)/365, number.people = 20, iterations = 10, ...){
  # ... allows 'k' to be passed as an argument to the statistics function
  # the aim is to do an MC simulation for an event tested with statistics function
  # and return a vector of 1's and 0's where 1 - event happened in a given room
  # 0 - event hasn't happened
  # if probabilities do not sum up to 1 (with a possibility
  # of numerical instability) return error
  if (abs(sum(birthday.probs)-1)>1e-20){
    stop('Probabilities do not add up to 1')
  }
  # declare list, where for each room there will be an indicator
  # telling whether a particular event tested with "statistics" fuction
  # has happened in a given room
  room.statistics <- list()
  for (i in 1:iterations){
    # create sample room and populate it with people
    # the vector contains days of birth of each person in a room
    room <- sample(1:length(birthday.probs), size = number.people, replace = T, prob = birthday.probs)
    # test whether a particular event has happened in a room
    # and store it in a list
    room.statistics[[i]]<-statistics(room, ...)
  }
  return(room.statistics)
}

simulate.room.d<- function(statistics, d = 365, number.people = 20, iterations = 10, ...){
  # ... allows 'k' to be passed as an argument to the statistics function
  # uniform birthday distribution
  # the aim is to do an MC simulation for an event tested with statistics function
  # and return a vector of 1's and 0's where 1 - event happened in a given room
  # 0 - event hasn't happened

  # declare list, where for each room there will be an indicator
  # telling whether a particular event tested with "statistics" fuction
  # has happened in a given room
  room.statistics <- list()
  for (i in 1:iterations){
    # create sample room and populate it with people
    # the vector contains days of birth of each person in a room
    room <- sample(1:d, size = number.people, replace = T)
    # test whether a particular event has happened in a room
    # and store it in a list
    room.statistics[[i]]<-statistics(room, ...)
  }
  return(room.statistics)
}

bdays.summary <- function(x){
  # helper function to change the format from a vector with numbers
  # to a summary table with dates of birth listed in the first column
  # and number of people in a room who are born at that day in the 2nd
  as.data.frame(unclass(rle(sort(x))))[,2:1]
}

same.day.bday <- function(x, k = 2){
  # test function
  # test whether there is at least one pair of people with the same bday
  # in a room
  x.summary <- bdays.summary(x)
  if (any(x.summary[,2]>=k))
    return(1)
  else
    return(0)
}

no.same.day.bday <- function(x, k = 2){
  # test function
  # test whether no two people share the same bday in a room
  return(1-same.day.bday(x, k = k))
}


summary.bday.simulation<- function(statistics, birthday.probs = rep(1,365)/365, number.people = 20, iterations = 10, boot.iter = 1e3, ...){
  # ... allows 'k' to be passed as an argument to the statistics function
  # main function - summary
  # returns 95% confidence interval from bootstrap
  mean.special<- function(data, indices){
    # helper function to pass "mean" to boot
    d <- data[indices]
    return(mean(d))
  }
  # do MC simulation of the event tested with statistics function
  many.rooms <- simulate.room(statistics, birthday.probs = birthday.probs, number.people = number.people, iterations = iterations, ...)
  # bootstrap the result
  result <- boot(data = unlist(many.rooms), statistic = mean.special, R = boot.iter)
  # return 95% confidence interval
  return(boot.ci(result, conf = 0.95, type = "basic"))
}

summary.bday.simulation.d<- function(statistics, d = 365, number.people = 20, iterations = 10, boot.iter = 1e3, ...){
  # ... allows 'k' to be passed as an argument to the statistics function
  # uniform birthday distribution
  # main function - summary
  # returns 95% confidence interval from bootstrap
  mean.special<- function(data, indices){
    # helper function to pass "mean" to boot
    d <- data[indices]
    return(mean(d))
  }
  # do MC simulation of the event tested with statistics function
  many.rooms <- simulate.room.d(statistics, d = d, number.people = number.people, iterations = iterations, ...)
  # bootstrap the result
  result <- boot(data = unlist(many.rooms), statistic = mean.special, R = boot.iter)
  # return 95% confidence interval
  return(boot.ci(result, conf = 0.95, type = "basic"))
}


##################
# From this point onwards we assume that probabilities of births are uniform
# throughout the year
##################

b.1 <- function(number.people = 20, num.same.day.bday = 2, num.days = 365){
  # calculate b_1 from the paper by ARRATIA, GOLDSTEIN and GORDON
  b1 <- choose(number.people, num.same.day.bday) *
    (choose(number.people, num.same.day.bday) - choose(number.people-num.same.day.bday, num.same.day.bday))*
    num.days^(2-2*num.same.day.bday)
  return(b1)
}

b.2 <- function(number.people = 20, num.same.day.bday = 2, num.days = 365){
  # calculate b_2 from the paper by ARRATIA, GOLDSTEIN and GORDON
  a <- choose(number.people, num.same.day.bday)
  total = 0
  for (j in seq.alt(1,num.same.day.bday-1)){
    c <- choose(number.people-num.same.day.bday, num.same.day.bday-j)
    e <- choose(num.same.day.bday, j)
    p <- num.days^(1+j-2*num.same.day.bday)
    total = total + c*e*p*a
  }
  return(total)
}

find.lambda <- function(number.people = 20, num.same.day.bday = 2, num.days = 365){
  # calculate lambda=E[W]
  lambda <- choose(number.people, num.same.day.bday) * num.days^(1-num.same.day.bday)
  return(lambda)
}

tv.bound.no.occurence <- function(number.people = 20, num.same.day.bday = 2, num.days = 365){
  # calculate the 2nd bound from theorem 1:
  # |P(W=0)-exp(-lambda)|<(b_1+b_2)exp(-lambda)/lambda
  # if one wants to extend this to a general bound on Total Variation
  # then use the first bound from theorem 1, i.e. multiply the result
  # of this function by 2
  lambda <- find.lambda(number.people, num.same.day.bday, num.days)
  b1 <- b.1(number.people, num.same.day.bday, num.days)
  b2 <- b.2(number.people, num.same.day.bday, num.days)
  bound <- (b1+b2)*(1-exp(-lambda))/lambda
  return(bound)
}
