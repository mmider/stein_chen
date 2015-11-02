source("plotting_functions.R")

n1 <- seq(11,31,by=2)
graph.n(n1)

n2 <- seq(40, 60, by=2)
graph.n(n2, k=3)

d3 <- seq(150,750,by=50)
graph.d(d=d3)

n3 <- seq(100, 5000, 100)
graph.nd(n1)
graph.nd(n3)

for (i in 1:length(n3)) {
  print(n3[i])
  d <- d.n.const.ratio(n3[i], 2)
  print(Bday_MC(c(0,0,0), n3[i], 2, d))
}  

graph.nd.no.sim(n3)
