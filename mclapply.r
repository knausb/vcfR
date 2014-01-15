

library(parallel)

setseed(1)
x <- matrix(as.integer(runif(20, min=1, max=20)), ncol=2)

as.vector(x)
as.list(x)

mclapply(x, sum)




# EOF.
