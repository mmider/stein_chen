######## THE LENGTH OF THE LONGEST HEAD RUN ######## 

######## FUNCTIONS:

require(boot)

# Indicator function for the presence of a run of t-heads in 
# sequence of coin tosses.

num.clusters.exceeding.t = function(coin.tosses, t, comp.pois = F){
  cluster.sizes = rle(coin.tosses)
  
  # We check if the first realisation of our sequence of tosses is 
  # head or tail. 
  
  if (cluster.sizes$values[1] == 0){
    # Selecting the even positions
    even.pos = (1:(length(cluster.sizes$values)%/% 2)) * 2
    head.cluster.sizes <- cluster.sizes$lengths[even.pos]
  }
  
  else{
    # Selecting the odd positions
    odd.pos <- (1:(ceiling(length(cluster.sizes$values) / 2 ))) * 2 - 1
    head.cluster.sizes <- cluster.sizes$lengths[odd.pos]
  }
  
  # Not considering clumbs (comp.pois == F)
  if (comp.pois == F){
    W <- sum(head.cluster.sizes >= t)
    return(W)
  }
  
  # Considering clumbs (comp.pois == T)
  else{
    U <- sum(head.cluster.sizes[head.cluster.sizes>=t]-t+1)
    return(U)
  }
}


# Auxiliary function to have a proper indicator function
exist.run.exceeding.t = function(vec){
  return(1*(vec > 0))
}


# Function to perform MC simulation for a given statistics

head.run.MC = function(test.function, tosses.length = 2047, t = 14, p = c(0.5, 0.5),
                       iterations = 1e2, comp.pois = F, boot.iter = 1e2) {
  
  
  #Auxiliary function to pass "mean" to boot
  mean.special = function(data, indices){
    d <- data[indices]
    return(mean(d))
  }
  
  counts = rep(0, iterations)
  size = tosses.length + t
  
  for (i in 1:iterations){
    coin.tosses = sample(c(0,1), size = size, replace = T, prob = p)
    if (comp.pois == F)
      var = num.clusters.exceeding.t(coin.tosses, t)
    else
      var = num.clusters.exceeding.t(coin.tosses, t, T)
    counts[i] = test.function(var)
  }
  
  result = boot(data = counts, statistic = mean.special, R = boot.iter)
  list( CI = boot.ci(result, conf = 0.95, type = "basic"), estimate = mean(counts) )
}


# Function to evluate lambda = E(W)
find.lambda = function(n, p, t) {
  (p ^ t)*((n - 1)*( 1 - p) + 1)
}


# Function evaluate (almost) b1:

b.1 = function(n, p, t, lambda) {
  b1 = ((lambda^2)*(2*t + 1)/n) + 2*lambda*p^t
  return(b1)
}


# Function to evaluate Stein-Chen boundary

tv.bound.no.occurence = function(n, p, t) {
  lambda = find.lambda(n, p, t)
  prob = exp(-lambda)
  b1 = b.1(n, p, t, lambda)
  bound = b1*(1 - exp(-lambda))/lambda
  bounds = abs(1 - c(prob - bound, prob + bound))
  return(bounds)
}
