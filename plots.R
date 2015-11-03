source("plotting_functions.R")

n1 <- seq(11,31,by=2)
pn.n1 <- graph.n(n1)
pn.n1$plot

n2 <- seq(40, 60, by=2)
pn.n2 <- graph.n(n2, k=3)
pn.n2$plot

d3 <- seq(150,750,by=50)
pd.d3 <- graph.d(d=d3)
pd.d3$plot

d4 <- seq(400,800, by=50)
pd.d4 <- graph.d(n = 30, k=2, d=d4)
pd.d4$plot

d4 <- seq(400,800, by=50)
pd.50.3.d4 <- graph.d(n = 50, k=3, d=d4)
pd.50.3.d4$plot


n3 <- seq(100, 5000, 100)
pnd.n1 <- graph.nd(n1)
pnd.n1$plot

for (i in 1:length(n3)) {
  print(n3[i])
  d <- d.n.const.ratio(n3[i], 2)
  print(Bday_MC(c(0,0,0), n3[i], 2, d))
}  

graph_list <- list(pn.n1, pn.n2, pnd.n1)
graph_list_plot <- list(pn.n1$plot, pn.n2$plot, pn.n1$plot)

nGraphsCol(graph_list_plot)


graph.nd.no.sim(n3)

library(MASS)
k <- 2
cbind(hills, k = k)
a <- c(2,3,4)
b <- c("a", "b", "c")

data.frame(a, b, k)

# graph.n, graph.d, graph.nd with faceting
pn.n1 <- graph.n(n1)
pn.n1$plot

pn.n2 <- graph.n(n2, k=3)
pn.n2$plot

graph_list <- list(pn.n1, pn.n2)
graph_list

graph_list_d <- list(pd.d3, pd.d4, pd.50.3.d4)

graph.n_facet(graph_list, "k")

graph.d_facet(graph_list_d, "n")



