library("ggplot2")
library("coda")
library("rjags")

# Current working directory of R
getwd()

# Change working directory of R
setwd("/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/R code/")

x <- read.csv("hw5.csv")
x <- x[, 2:4]

n <- ncol(x)
m <- nrow(x)

k <- 3

jags_model_string <- "model {
    for (i in 1:n*m) {
        x[i] ~ dpois(theta[group[i]])
    }  
    
    for (j in 1:n) {
        theta[j] ~ dgamma(k, invsigma)
    }

    invsigma ~ dexp(1.0)
}"

x_jags <- array(dim = c(m * n, 2))

colnames(x_jags) <- c("x", "group")

x_jags[1:100, 1] <- x[, 1]
x_jags[1:100, 2] <- 1.0
x_jags[101:200, 1] <- x[, 2]
x_jags[101:200, 2] <- 2.0
x_jags[201:300, 1] <- x[, 3]
x_jags[201:300, 2] <- 3.0

data_jags <- list(
    x = x_jags[, 1],
    group = x_jags[, 2],
    n = n,
    m = m,
    k = k
)

model <- jags.model(
    textConnection(jags_model_string),
    data = data_jags,
    n.chains = 3,
)

update(model, 100) # burn -in period

params <- c("theta", "invsigma")

model_simulation <- coda.samples(
    model = model,
    variable.names = params,
    n.iter = 10000
)

summary(model_simulation)