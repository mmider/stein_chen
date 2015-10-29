####### First Case: SIMPLE BIRTHDAY PROBLEM.
####### UNIFORM BIRTHDAYS, 2-WAY COINCIDENCE

# n individuals
# birthday of individuals are indepedendent 
# d days


# Function to evaluate the probability of no coincidence
prob.no.coincidence = function(n, d) {
  temp1 = rep(1, n - 1)
  temp2 = seq(from = 1, to = n - 1)/d
  tempseq = temp1 - temp2
  return(prod(tempseq))
}


####### Second Case: BIRTHDAY PROBLEM.
####### UNIFORM BIRTHDAYS, K-WAY COINCIDENCE

