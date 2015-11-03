library(Rcpp)
library(ggplot2)
library(grid)
sourceCpp("Bday_prob.cpp")
source("Bday_problem.R")

# helper function to calculate d as a function of n, keeping the ratio n^(k+1)/d^k constant 
# (constant set by ratio argument)
d.n.const.ratio <- function(n, k, ratio = 0.1) {
  ((n^(1+k))/ratio)^(1/k)
}

# plots theoretical Stein-Chen limits and simulation CIs, for probability that no k people share a birthday,
# for a range of n (number of people), given a fixed d (days in the year)
graph.n <- function(n, k = 2, d = 365, ...) {
  # n is a vector to loop over - each element is the number of people in the room for that iteration
  # k is the number of people who share a birthday
  # d is the number of days from which the birthdays can be picked
  # assume uniform distribution of birthdays
  
  # initialise vectors for storing bound and CI values
  len_n <- length(n)
  upper_CI <- numeric(len_n)
  lower_CI <- numeric(len_n)
  upper_lim <- numeric(len_n)
  lower_lim <- numeric(len_n)
  bound <- c(rep("upper", len_n), rep("lower", len_n))
  
  for (i in 1:len_n) {
    # bootstrap simulation of probability that no k of the n[i] people share a birthday
    # and store upper and lower 95% confidence interval bounds
    suma <- Bday_MC(c(0,0,0), n[i], k, d, ...)
    upper_CI[i] <- suma[3]
    lower_CI[i] <- suma[1]
    
    # calculate Poisson approximation that no k of the n[i] people share a birthday
    lambda <- find.lambda(n[i], k, d)
    w.zero <- exp(-lambda)
    
    # calculate and store Stein-Chen bounds for the Poisson approximation
    err <- tv.bound.no.occurence(n[i], k, d)
    upper_lim[i] <- min(1, w.zero + err)
    lower_lim[i] <- max(0, w.zero - err)
  }
  
  # to determine if Stein-Chen bounds are 'good', want to check they are no more than 0.05 away from simulated probability
  upper_CI_0.05 <- min(1, upper_CI + 0.05)
  lower_CI_0.05 <- max(0, lower_CI - 0.05)
  
  # put calculations in data frames, in suitable format for passing to ggplot2
  df_CI <- data.frame(n=rep(n,2), y=c(upper_CI, lower_CI), d=d, k=k)
  df_lim <- data.frame(n=rep(n,2), y=c(upper_lim, lower_lim), bound=bound, d=d, k=k)
  df_0.05 <- data.frame(n=rep(n,2), y=c(upper_CI_0.05, lower_CI_0.05), bound=bound, d=d, k=k)
  
  # create and return plot
  q <- ggplot(NULL, aes(n, y)) + geom_line(data=df_lim, aes(group=bound, colour="bound")) + geom_point(data=df_lim, aes(colour="bound")) 
  q <- q + geom_point(data=df_CI, aes(colour="CI")) + geom_line(data=df_CI, aes(group=n, colour="CI"))
  q <- q + geom_line(data=df_0.05, aes(group=bound, colour="CI +/- 0.05"), linetype="dashed") + geom_point(data=df_0.05, aes(colour = "CI +/- 0.05"))
  return(list(plot = q, df_CI = df_CI, df_lim = df_lim, df_0.05 = df_0.05))
}

# plots theoretical Stein-Chen limits and simulation CIs, for probability that no k of n people share a birthday,
# for a range of d (number of days from which birthdays can be picked)
graph.d <- function(n = 23, k = 2, d, ...) {
  # n is the number of people in the room 
  # k is the number of people who share a birthday
  # d is a vector to be looped over, with each element being 
    # the number of days from which the birthdays can be picked in that iteration
  # assume uniform distribution of birthdays
  
  # initialise vectors for storing bound and CI values
  len_d <- length(d)
  upper_CI <- numeric(len_d)
  lower_CI <- numeric(len_d)
  upper_lim <- numeric(len_d)
  lower_lim <- numeric(len_d)
  bound <- c(rep("upper", len_d), rep("lower", len_d))
  
  for (i in 1:len_d) {
    # bootstrap simulation of probability that no k of the n people share a birthday, from d[i] possible birthdays
    # and store upper and lower 95% confidence interval bounds
    suma <- Bday_MC(c(0,0,0), n, k, d[i], ...)
    upper_CI[i] <- suma[3]
    lower_CI[i] <- suma[1]

    # calculate Poisson approximation that no k of the n people share a birthday
    lambda <- find.lambda(n, k, d[i])
    w.zero <- exp(-lambda)
    
    # calculate and store Stein-Chen bounds for the Poisson approximation
    err <- tv.bound.no.occurence(n, k, d[i])
    upper_lim[i] <- min(1, w.zero + err)
    lower_lim[i] <- max(0, w.zero - err)
    }

  # to determine if Stein-Chen bounds are 'good', want to check they are no more than 0.05 away from simulated probability
  upper_CI_0.05 <- min(1, upper_CI + 0.05)
  lower_CI_0.05 <- max(0, lower_CI - 0.05)
  
  # put calculations in data frames, in suitable format for passing to ggplot2
  df_CI <- data.frame(d=rep(d,2), y=c(upper_CI, lower_CI), n=n, k=k)
  df_lim <- data.frame(d=rep(d,2), y=c(upper_lim, lower_lim), bound=bound, n=n, k=k)
  df_0.05 <- data.frame(d=rep(d,2), y=c(upper_CI_0.05, lower_CI_0.05), bound=bound, n=n, k=k)
  
  # create and return plot
  q <- ggplot(NULL, aes(d, y)) + geom_line(data=df_lim, aes(group=bound, colour="bound")) + geom_point(data=df_lim, aes(colour="bound")) 
  q <- q + geom_point(data=df_CI, aes(colour="CI")) + geom_line(data=df_CI, aes(group=d, colour="CI"))
  q <- q + geom_line(data=df_0.05, aes(group=bound, colour="CI +/- 0.05"), linetype="dashed") + geom_point(data=df_0.05, aes(colour = "CI +/- 0.05"))
  return(list(plot = q, df_CI = df_CI, df_lim = df_lim, df_0.05 = df_0.05))
}

# plots theoretical Stein-Chen limits and simulation CIs, for probability that no k people share a birthday,
# for a range of n (number of people), with d (days in the year) s.t. the ratio n^(k+1)/d^k is kept constant,
# (the value of that constant is the ratio argument)
graph.nd <- function(n, k = 2, ratio = 0.1, ...) {
  # n is a vector to loop over - each element is the number of people in the room for that iteration
  # k is the number of people who share a birthday
  # ratio is the value that n^(k+1)/d^k is to be kept constant at
  # assume uniform distribution of birthdays
  
  # initialise vectors for storing bound and CI values
  len_n <- length(n)
  
  upper_CI <- numeric(len_n)
  lower_CI <- numeric(len_n)
  upper_lim <- numeric(len_n)
  lower_lim <- numeric(len_n)
  bound <- c(rep("upper", len_n), rep("lower", len_n))
  
  for (i in 1:len_n) {
    # find value of d which keeps the ratio n^(k+1)/d^k constant (that constant is the ratio argument)
    d <- d.n.const.ratio(n[i], k, ratio)
    
    # bootstrap simulation of probability that no k of the n[i] people share a birthday, from d possible birthdays
    # and store upper and lower 95% confidence interval bounds
    suma <- Bday_MC(c(0,0,0), n[i], k, d, ...)
    upper_CI[i] <- suma[3]
    lower_CI[i] <- suma[1]
    
    # calculate Poisson approximation that no k of the n people share a birthday
    lambda <- find.lambda(n[i], k, d)
    w.zero <- exp(-lambda)
    
    # calculate and store Stein-Chen bounds for the Poisson approximation
    err <- tv.bound.no.occurence(n[i], k, d)
    upper_lim[i] <- min(1, w.zero + err)
    lower_lim[i] <- max(0, w.zero - err)
  }
  
  # to determine if Stein-Chen bounds are 'good', want to check they are no more than 0.05 away from simulated probability
  upper_CI_0.05 <- min(1, upper_CI + 0.05)
  lower_CI_0.05 <- max(0, lower_CI - 0.05)
  
  # put calculations in data frames, in suitable format for passing to ggplot2
  df_CI <- data.frame(n=rep(n,2), y=c(upper_CI, lower_CI), k=k)
  df_lim <- data.frame(n=rep(n,2), y=c(upper_lim, lower_lim), bound=bound, k=k)
  df_0.05 <- data.frame(n=rep(n,2), y=c(upper_CI_0.05, lower_CI_0.05), bound=bound, k=k)
  
  # create and return plot
  q <- ggplot(NULL, aes(n, y)) + geom_line(data=df_lim, aes(group=bound, colour="bound")) + geom_point(data=df_lim, aes(colour="bound")) 
  q <- q + geom_point(data=df_CI, aes(colour="CI")) + geom_line(data=df_CI, aes(group=n, colour="CI"))
  q <- q + geom_line(data=df_0.05, aes(group=bound, colour="CI +/- 0.05"), linetype="dashed") + geom_point(data=df_0.05, aes(colour = "CI +/- 0.05"))
  return(list(plot = q, df_CI = df_CI, df_lim = df_lim, df_0.05 = df_0.05))
}

# like graph.nd, but without the simulation (since this takes a long time to run for large n)
graph.nd.no.sim <- function(n, k = 2, ratio = 0.1, ...) {
  # n is a vector to loop over - each element is the number of people in the room for that iteration
  # k is the number of people who share a birthday
  # ratio is the value that n^(k+1)/d^k is to be kept constant at
  # assume uniform distribution of birthdays
  
  # initialise vectors for storing bound and CI values
  len_n <- length(n)
#  upper_CI <- numeric(len_n)
#  lower_CI <- numeric(len_n)
  upper_lim <- numeric(len_n)
  lower_lim <- numeric(len_n)
  bound <- c(rep("upper", len_n), rep("lower", len_n))
  
  for (i in 1:len_n) {
    # find value of d which keeps the ratio n^(k+1)/d^k constant (that constant is the ratio argument)
    d <- d.n.const.ratio(n[i], k, ratio)
    
    # bootstrap simulation of probability that no k of the n[i] people share a birthday, from d possible birthdays
    # and store upper and lower 95% confidence interval bounds
#    suma <- Bday_MC(c(0,0,0), n[i], k, d, ...)
#    upper_CI[i] <- suma[3]
#    lower_CI[i] <- suma[1]
    
    # calculate Poisson approximation that no k of the n people share a birthday
    lambda <- find.lambda(n[i], k, d)
    w.zero <- exp(-lambda)
    
    # calculate and store Stein-Chen bounds for the Poisson approximation
    err <- tv.bound.no.occurence(n[i], k, d)
    upper_lim[i] <- min(1, w.zero + err)
    lower_lim[i] <- max(0, w.zero - err)
  }
  
  # to determine if Stein-Chen bounds are 'good', want to check they are no more than 0.05 away from simulated probability
#  upper_CI_0.05 <- upper_CI + 0.05
#  lower_CI_0.05 <- lower_CI - 0.05
  
  # put calculations in data frames, in suitable format for passing to ggplot2
#  df_CI <- data.frame(n=rep(n,2), y=c(upper_CI, lower_CI))
  df_lim <- data.frame(n=rep(n,2), y=c(upper_lim, lower_lim), bound=bound)
#  df_0.05 <- data.frame(n=rep(n,2), y=c(upper_CI_0.05, lower_CI_0.05), bound=bound)
  
  # create and return plot
  q <- ggplot(NULL, aes(n, y)) + geom_line(data=df_lim, aes(group=bound, colour="bound")) + geom_point(data=df_lim, aes(colour="bound")) 
#  q <- q + geom_point(data=df_CI, aes(colour="CI")) + geom_line(data=df_CI, aes(group=n, colour="CI"))
#  q <- q + geom_line(data=df_0.05, aes(group=bound, colour="CI +/- 0.05"), linetype="dashed") + geom_point(data=df_0.05, aes(colour = "CI +/- 0.05"))
  q
}

# function to print n graphs in a column (generic)
nGraphsCol <- function(graphs) {
  # graphs is a list of graphs
  n <- length(graphs)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(n,1)))
  vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
  for (i in 1:n) {
    print(graphs[[i]], vp=vplayout(i,1))
  }
}

# function to print n graphs in a row (generic)
nGraphsRow <- function(graphs) {
  # graphs is a list of graphs
  n <- length(graphs)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,n)))
  vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
  for (i in 1:n) {
    print(graphs[[i]], vp=vplayout(1,i))
  }
}

# facets a list of graphs by "k" or "d"
# NB "k" or "d" must be different for each graph
graph.n_facet <- function(graph_list, facet_by) {
  # each item in graph_list must be the output from one of graph.n or graph.nd
  # facet_by must be "k" or "d"
  
  len <- length(graph_list)
  
  df_lim <- graph_list[[1]]$df_lim
  df_CI <- graph_list[[1]]$df_CI
  df_0.05 <- graph_list[[1]]$df_0.05
  
  for (i in 2:len) {
    df_lim <- rbind(df_lim, graph_list[[i]]$df_lim)
    df_CI <- rbind(df_CI, graph_list[[i]]$df_CI)
    df_0.05 <- rbind(df_0.05, graph_list[[i]]$df_0.05)
  }
  
  q <- ggplot(NULL, aes(n, y)) + geom_line(data=df_lim, aes(group=bound, colour="bound")) + geom_point(data=df_lim, aes(colour="bound")) 
  q <- q + geom_point(data=df_CI, aes(colour="CI")) + geom_line(data=df_CI, aes(group=n, colour="CI"))
  q <- q + geom_line(data=df_0.05, aes(group=bound, colour="CI +/- 0.05"), linetype="dashed") + geom_point(data=df_0.05, aes(colour = "CI +/- 0.05"))
  
  if (facet_by == "k") q <- q + facet_grid(.~k, scales = "free")
  if (facet_by == "d") q <- q + facet_grid(.~d, scales = "free")
  q
}

# facets a list of graphs by "k" or "n"
# NB "k" or "n" must be different for each graph
graph.d_facet <- function(graph_list, facet_by) {
  # each item in graph_list must be the output from one of graph.n or graph.nd
  # facet_by must be "k" or "n"
  
  len <- length(graph_list)
  
  df_lim <- graph_list[[1]]$df_lim
  df_CI <- graph_list[[1]]$df_CI
  df_0.05 <- graph_list[[1]]$df_0.05
  
  for (i in 2:len) {
    df_lim <- rbind(df_lim, graph_list[[i]]$df_lim)
    df_CI <- rbind(df_CI, graph_list[[i]]$df_CI)
    df_0.05 <- rbind(df_0.05, graph_list[[i]]$df_0.05)
  }
  
  q <- ggplot(NULL, aes(d, y)) + geom_line(data=df_lim, aes(group=bound, colour="bound")) + geom_point(data=df_lim, aes(colour="bound")) 
  q <- q + geom_point(data=df_CI, aes(colour="CI")) + geom_line(data=df_CI, aes(group=d, colour="CI"))
  q <- q + geom_line(data=df_0.05, aes(group=bound, colour="CI +/- 0.05"), linetype="dashed") + geom_point(data=df_0.05, aes(colour = "CI +/- 0.05"))
  
  if (facet_by == "k") q <- q + facet_grid(.~k, scales = "free")
  if (facet_by == "n") q <- q + facet_grid(.~n, scales = "free")
  q
}
