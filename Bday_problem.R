# comment your code!
require('boot')
simulate.room<- function(statistics, birthday.probs = rep(1,365)/365, number.people = 20, iterations = 10){
  if (abs(sum(birthday.probs)-1)>1e-20){
    stop('Probabilities do not add up to 1')
  }
  room.statistics <- list()
  for (i in 1:iterations){
    room <- sample(1:length(birthday.probs), size = number.people, replace = T, prob = birthday.probs)
    room.statistics[[i]]<-statistics(room)
  }
  return(room.statistics)
}

bdays.summary <- function(x){
  as.data.frame(unclass(rle(sort(x))))[,2:1]
}

same.day.bday <- function(x, k = 2){
  x.summary <- bdays.summary(x)
  if (any(x.summary[,2]>=k))
    return(1)
  else
    return(0)
}

no.same.day.bday <- function(x, k = 2){
  return(1-same.day.bday(x, k = k))
}


summary.bday.simulation<- function(statistics, birthday.probs = rep(1,365)/365, number.people = 20, iterations = 10, boot.iter = 1e3){
  mean.special<- function(data, indices){
    d <- data[indices]
    return(mean(d))
  }  
  many.rooms <- simulate.room(statistics, birthday.probs = birthday.probs, number.people = number.people, iterations = iterations)
  result <- boot(data = unlist(many.rooms), statistic = mean.special, R = boot.iter)
  return(boot.ci(result, conf = 0.95, type = "basic"))
}
