source('Head_run.R')


tosses.length = 2047
t = 14
iterations = 1e3
boot.iter = 1e2
counts <- rep(0, iterations)
t = 14

# using W
head.run.MC(exist.run.exceeding.t, tosses.length, t, iterations, comp.pois = F, boot.iter)

# using U
head.run.MC(exist.run.exceeding.t, tosses.length, t, iterations, comp.pois = T, boot.iter)


