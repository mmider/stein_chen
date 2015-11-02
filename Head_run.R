require(boot)

num.clusters.exceeding.t <- function(coin.tosses, t, comp.pois = F){
  cluster.sizes <- rle(coin.tosses)
  if (cluster.sizes$values[1] == 0){
    # pick out all the even positions
    even.pos <- (1:(length(cluster.sizes$values)%/% 2)) * 2
    head.cluster.sizes <- cluster.sizes$lengths[even.pos]
  }
  else{
    # pick out odd positions
    odd.pos <- (1:(ceiling(length(cluster.sizes$values) / 2 ))) * 2 - 1
    head.cluster.sizes <- cluster.sizes$lengths[odd.pos]
  }
  if (comp.pois == F){
    W <- sum(head.cluster.sizes >= t)
    return(W)
  }
  else{
    U <- sum(head.cluster.sizes[head.cluster.sizes>=t]-t+1)
    return(U)
  }
}

exist.run.exceeding.t<-function(vec){
  return(1*(vec > 0))
}



head.run.MC <- function(test.function, tosses.length = 2047, t = 14, iterations = 1e2, comp.pois = F, boot.iter = 1e2){
  mean.special<- function(data, indices){
    # helper function to pass "mean" to boot
    d <- data[indices]
    return(mean(d))
  }
  size = tosses.length + t
  for (i in 1:iterations){
    coin.tosses <- sample(c(0,1), size = size, replace = T)
    if (comp.pois == F)
      var <- num.clusters.exceeding.t(coin.tosses,t)
    else
      var <- num.clusters.exceeding.t(coin.tosses,t,T)
    counts[i]<-test.function(var)

  }
  result <- boot(data = counts, statistic = mean.special, R = boot.iter)
  list( CI = boot.ci(result, conf = 0.95, type = "basic"), estimate = mean(counts) )
}



