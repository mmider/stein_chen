# IMPORTANT
# set working directory to the folder containing Bday_problem.R:
source('Biological_example.R')

data <- import.DNA(512)
DNA <- data$data
mu <- data$mu
p <- data$p
num.pairs <- 200
window.size <- 21
string.length <-512
a <- 0.8069809
s <- 16
iterations <- 1e4
init.value <- 1

calculate.bound.thm3(string.length,string.length,window.size,s,a,mu,p, iterations,init.value)
calculate.bound.thm3(1024,1024,37,27,0.9,mu,p,1e4,1)

synthetic.DNA <- simulate.dna(mu,nrow(DNA),ncol(DNA))

# Main plots for report
plot.hang.hist(DNA,num.pairs,window.size,p,"Liverwort sequence -")
plot.hang.hist(synthetic.DNA,num.pairs,window.size,p, "Simulated DNA sequence -")
# should return 21:
# calc.t(a,string.length,p)