tvar2sim=function (sd = 1) 
{
    n <- 512
    ar1vec <- seq(from = -1.1, to = 0.5, length = n)
    ar2vec <- seq(from = -1.1, to = 0.5, length = n)
    x <- c(rnorm(2, mean = 0, sd = sd), rep(0, n - 2))
    for (i in 3:n) x[i] <- ar1vec[i] * x[i - 1] + ar2vec[i]*x[i-2] + rnorm(1, mean = 0, sd = sd)
    return(x)
}

