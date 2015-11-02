require(rootSolve)
require(Rcpp)
require(ggplot2)

sourceCpp('Biological_example.cpp')

import.DNA <- function(stripe.length = 512){
  # import data from file containing DNA of liverwort.
  # Data downloaded from GenBank: http://www.ncbi.nlm.nih.gov/nuccore/11640?report=fasta
  file.name <- 'MarchantiaPolymorpha.txt'
  # clean up the data
  data <- readChar(file.name, file.info(file.name)$size)
  data <- gsub('\n','',data)
  data <- gsub('\r','',data)
  # should be equal to 121024
  if (nchar(data) != 121024){
    stop('invalid number of bases, check the source file')
  }
  
  # convert to vector with entries - letters of DNA sequence
  data <- as.factor(strsplit(data, split="")[[1]])
  # prob. of seeing given letters
  mu <- summary(data)/length(data)
  # prob. of a match between any two letters
  p <- sum(mu^2)
  # number of stripes which form by splitting DNA strand into stipes of
  # size stripe.length each. The remaining letters are discarded
  num.stripes <- floor(length(data)/stripe.length)
  # return as a matrix, where in each column we have a stripe of DNA
  data <- matrix(data[1:(num.stripes * stripe.length)], ncol = num.stripes)
  return( list(data = data, mu = mu, p = p))
}

find.longest.runs<-function(DNA, num.pairs = 200, window.size = 21, prohibit.repeats = T){
  # We make num.pairs comparisons. Each time we draw at random two pairs
  # of our DNA stripes. The we slide the window on each of them and count
  # number of matches of letters. We do it for all possible combinations
  # of positions of two windows on two stripes and we store the largest number
  # of matches that we have found in such experiment. We repeat this procedure
  # num.pairs times yielding num.pairs maximal matches, which we then return

  # sample num.pairs different pairings of DNA
  sample.pairs.A <- sample(1:ncol(DNA), size = num.pairs, replace = F)
  sample.pairs.B <- sample(1:ncol(DNA), size = num.pairs, replace = F)
  # we can specify that we never want to compare the same pair twice:
  if (prohibit.repeats){
    while (any(sample.pairs.A-sample.pairs.B==0)){
      sample.pairs.A <- sample(1:ncol(DNA), size = num.pairs, replace = F)
      sample.pairs.B <- sample(1:ncol(DNA), size = num.pairs, replace = F)
    }
  }
  # store the sampled pairs of DNA stripes
  DNA.A <- DNA[,sample.pairs.A]
  DNA.B <- DNA[,sample.pairs.B]
  # declare vector, where maximal matches will be stored
  longest.runs <- rep(0, num.pairs)
  for (i in 1:num.pairs){
    # this is a Cpp function which does all possible comparisons and returns
    # the maximum number of matches
    longest.runs[i] <- longest_common_run(DNA.A[,i], DNA.B[,i], window.size)
  }
  return(longest.runs)
}

predict.freq.max.run <- function(observed.max, stripe.length, window.size, prob.match){
  # use chen-stein poisson approximation to esimate probability
  # that observed.max will be the maximum number of matches
  # we will see in our experiment - just follow equation (31) in paper
  s <- observed.max
  t <- window.size
  n <- stripe.length
  p <- prob.match
  estim.prob <- exp(-(n-t)^2*((s+1)/t-p)*dbinom(s+1,t,p)) - 
    exp(-(n-t)^2*(s/t-p)*dbinom(s,t,p))
  return(estim.prob)
}

adjust.empirical.probs <- function(empirical.probs, window.size){
  # technical adjustment - fill in vector with zeros
  empirical.adjusted <- c(empirical.probs, rep(0,(window.size - length(empirical.probs))))
  return(empirical.adjusted)
}

arcsine.tr <- function(x){
  return(asin(sqrt(x)))
}

histogram.centre <- function(n,t,p,lower.guess,upper.guess){
  # Functions to calculate centre of hanging histograms
  # (equivalently most frequent maximum run length):
  root.function <- function(s,n,t,p){
    # function from Erdos-Renyi paper
    out <- n^2 * (s/t-p)*gamma(1+t)/(gamma(1+s)*gamma(1+t-s))*p^s*(1-p)^(t-s)-1
    return(out)
  }
  root.fun.wrap <- function(s){
    root.function(s,n,t,p)
  }
  root <- floor(uniroot(root.fun.wrap, c(lower.guess,upper.guess))$root)
  return(root)
}

simulate.dna <- function(mu, string.len, total.strings){
  # simulates synthetic DNA data
  sim.DNA <- matrix(sample(c("A", "C", "G", "T"),
                           size = string.len * total.strings,
                           replace = T),
                    ncol = total.strings)
  return(sim.DNA)
}

plot.hang.hist <- function(DNA,num.pairs, window.size, p, title.add = "", prohibit.repeats = T){
  # draws hanging histograms
  n <- nrow(DNA)
  longest.runs <- find.longest.runs(DNA,num.pairs, window.size, prohibit.repeats)
  empirical.probs <- adjust.empirical.probs(tabulate(longest.runs)/num.pairs, window.size)
  
  hist.centre <- histogram.centre(n,window.size,p,window.size/2,window.size - 1)
  
  x.range <- (hist.centre - 3):(hist.centre+3)
  n.temp <- length(x.range)
  obs.freq <- empirical.probs[x.range]
  pred.freq <- predict.freq.max.run(x.range,n,window.size,p)
  
  df <- data.frame(x = rep(x.range,2),
                   y = c(arcsine.tr(pred.freq),arcsine.tr(pred.freq) - arcsine.tr(obs.freq)),
                   legend = c(rep("predicted", n.temp),rep("predicted - observed", n.temp))
  )
  df2 <- data.frame(
    x = x.range,
    y.pred = arcsine.tr(pred.freq),
    y.obs = arcsine.tr(obs.freq)
  )

  df3 <- data.frame(
    x = c(x.range-0.5, x.range[length(x.range)]+0.5),
    line.up = rep(1/(2*sqrt(num.pairs))*2.56, length(x.range)+1),
    line.mid = rep(0,length(x.range)+1),
    line.down = rep(-1/(2*sqrt(num.pairs))*2.56, length(x.range)+1)
  )
  df4 <- data.frame(
    x = rep(x.range[1]-0.55,2),
    y = c(df3$line.down[1],df3$line.up[1])
  )
  pt <- ggplot() +
    geom_rect(data = df2, aes(xmin = y.pred-y.obs, xmax = y.pred, ymin = x-0.5, ymax = x + 0.5, alpha = 0.1)) + 
    geom_linerange(data = df, aes(y, ymin = x - 0.5, ymax = x+0.5, colour = legend, size = 2))  +
    geom_line(data = df3, aes(x = line.up, y = x), linetype = 'dotted') +
    geom_line(data = df3, aes(x = line.mid, y = x), linetype = 'dotted') + 
    geom_line(data = df3, aes(x = line.down, y = x), linetype = 'dotted') + coord_flip() + guides(size = F, alpha = F) +
    geom_line(data = df4, aes(x = y, y = x), linetype = 'longdash', colour = 'cadetblue')+
    annotate("text", x = 0, y = x.range[1]-1, label = "+/-2.57s.e.")+
    labs(title = paste(title.add, "maximum matches between",window.size,"segments in", num.pairs, "pairs of ", n, "letter sequences"),
         y = paste("Maximum number of matches - out of possible", window.size),
         x = expression(paste("arcsin(", sqrt("freq"),")")))
  return(pt)
}

#######
# Calculate bound for Theorem 3 in paper about Erdos-Renyi Law in Distribution
# It seem to give a correct bound, but seems like a bit too much to explain
# in a report. Also one needs to add a bound from Theorem 1, but I haven't calculated
# that one.
#######
# helper function (ignore)
calc.t <- function(a,string.length,p){
  log((a-p)*string.length^2)/H(a,p)
}


H <- function(a,p){
  kullback.leibler.div <- a * log(a/p) + (1-a) * log((1-a)/(1-p))
  return(kullback.leibler.div)
}

G <- function(a, mu, p){
  c <- sum(mu^3)/p
  b.0 <- 1 - p*(1-c)^2/(2*a*(c-p))*((1+4*a*(1-a)*(c-p)/(p*(1-c)^2))^0.5-1)
  out <- a * H(b.0,c) + (1-a)*H((a-a*b.0)/(1-a),(p-c*p)/(1-p))
  return(out)
}

theta <- function(a, mu){
  p <- sum(mu^2)
  out <- G(a,mu,p) / H(a,p)
  return(out)
}

M.matrix <- function(mu, gamma){
  M <- sqrt(mu %o% mu) + diag(mu * (exp(gamma)-1))
  return(M)
}
f <- function(mu, gamma){
  lambda <- eigen(M.matrix(mu,gamma))$values[1]
  return(log(lambda))
}

increment.temp <- function(beta, alpha = 1.001, mu){
  return(beta * alpha)
}

I.MCMC <- function(a,p, mu, iterations, init = 1){
  I.maximise <- function(gamma){
    return(a*gamma - f(mu, gamma))
  }
  # finds the maximum of a function I(a) = sup_gamma (a * gamma - f(gamma))
  # current temperature
  X.path <- rep(0,iterations)
  X.path[1] <- init
  beta <- 1
  for (i in 1:(iterations-1)){
    beta <- increment.temp(beta)
    X.prop <- rnorm(1, mean = X.path[i], sd = 1)
    alpha <- exp(beta * (I.maximise(X.prop) - I.maximise(X.path[i])) )
    if (runif(1)<alpha){
      X.path[i+1] = X.prop
    }
    else{
      X.path[i+1] = X.path[i]
    }
  } 
  return(list(path = X.path, beta = beta, value = I.maximise(X.path[iterations])))
}

calculate.bound.thm3 <- function(m,n,t,s,a, mu,p,iterations,init = 1){
  I <- I.MCMC(a,p,mu,iterations,init)$value
  bound <- (4*t^2 + 5*t*(m+n) + 7*t + 1)*
    (1-pbinom(s-1,t,p))+ 
       4 * t^2*m*n*exp(-2*t*I)+
       2 * t * m * n^2 * exp(-t * (1 + theta(a, mu)) * H(a,p))+
       2 * t * n * m^2 * exp(-t * (1 + theta(a, mu)) * H(a,p))
  return(bound)
}


