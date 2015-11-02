######## THE LENGTH OF THE LONGEST HEAD RUN ######## 

######## Plotting boundaries and simulations.


source("HeadRuns_mod.R")



######## SIMULATION 1 ########



#### Setting parameters.

# In this case I fixed t = 8, and I'll let n grow

t = 8
p = 0.5
boot.iter = 1e3
iterations = 1e4

# Creating a sequence of n's.
seq.n = seq(from = 20, to = 40)

# Creating conf.intervals for the relative value
bounds = matrix(NA, nrow = length(seq.n), ncol = 2)
for(i in 1:length(seq.n)) {
  bounds[i, ] = tv.bound.no.occurence(n = seq.n[i], p = p, t = t)
}

# Creating vector of MC and relative confidence interval matrix
valMC = c()
confint.valMC = matrix(NA, nrow = length(seq.n), ncol = 2)

# Monte Carlo Simulation & populating vector and matrix
for(i in 1:length(seq.n)) {
  sample = head.run.MC(exist.run.exceeding.t, tosses.length = seq.n[i], t = 8, 
                       p = c(0.5,0.5), iterations, comp.pois = F, boot.iter)
  valMC[i] = sample$estimate
  confint.valMC[i, ] = sample$CI$basic[4:5]
}


# Plotting SC boundaries, MC simulations, and MC confidence intervals
plot(seq.n, bounds[, 1], type = "l", col = "red", lty = 2, lwd = 2,
     main = expression(paste("Stein-Chen boundaries for P[", R[n], "> 8]")),
     ylab = "Probability", xlab = "n")
lines(seq.n, bounds[, 2], type = "l", col = "red", lty = 2, lwd = 2)
lines(seq.n, valMC, type = "p", col = "green", pch = 19)
lines(seq.n, confint.valMC[, 1], type = "l", col = "green", lty = 2, lwd = 1)
lines(seq.n, confint.valMC[, 2], type = "l", col = "green", lty = 2, lwd = 1)


# Adding legend: Notice that with locator(1) you need to 
# select where to put the legend with your mouse.
legend(locator(1),  
       c("Stein-Chen boundaries","MC simulation", "MC 95 % CI"),
       lty=c(2,NA,2), pch=c(NA,19,NA), merge=FALSE, 
       col = c("red","green", "green"), lwd = c(2, NA, 1),
       cex = 0.8)






######### SIMULATION 2 ############


#### Setting parameters.

# In this case I fixed t = 8, and I'll let n grow

t = 4
p = 0.7
p_vec =
boot.iter = 1e3
iterations = 1e4

# Creating a sequence of n's.
seq.n = seq(from = 20, to = 40)

# Creating conf.intervals for the relative value
bounds = matrix(NA, nrow = length(seq.n), ncol = 2)
for(i in 1:length(seq.n)) {
  bounds[i, ] = tv.bound.no.occurence(n = seq.n[i], p = p, t = t)
}

# Creating vector of MC and relative confidence interval matrix
valMC = c()
confint.valMC = matrix(NA, nrow = length(seq.n), ncol = 2)

# Monte Carlo Simulation & populating vector and matrix
# (remember that when you pass the vector of probabilities 
# for the samplers, it in the position 1, the probability of 0)

for(i in 1:length(seq.n)) {
  sample = head.run.MC(exist.run.exceeding.t, tosses.length = seq.n[i], t = 8, 
                       p = c(0.3,0.7), iterations, comp.pois = F, boot.iter)
  valMC[i] = sample$estimate
  confint.valMC[i, ] = sample$CI$basic[4:5]
}


# Plotting SC boundaries, MC simulations, and MC confidence intervals
plot(seq.n, bounds[, 1], type = "l", col = "red", lty = 2, lwd = 2,
     main = expression(paste("Stein-Chen boundaries for P[", R[n], "> 8]")),
     ylab = "Probability", xlab = "n")
lines(seq.n, bounds[, 2], type = "l", col = "red", lty = 2, lwd = 2)
lines(seq.n, valMC, type = "p", col = "green", pch = 19)
lines(seq.n, confint.valMC[, 1], type = "l", col = "green", lty = 2, lwd = 1)
lines(seq.n, confint.valMC[, 2], type = "l", col = "green", lty = 2, lwd = 1)


# Adding legend: Notice that with locator(1) you need to 
# select where to put the legend with your mouse.
legend(locator(1),  
       c("Stein-Chen boundaries","MC simulation", "MC 95 % CI"),
       lty=c(2,NA,2), pch=c(NA,19,NA), merge=FALSE, 
       col = c("red","green", "green"), lwd = c(2, NA, 1),
       cex = 0.8)


  

       
