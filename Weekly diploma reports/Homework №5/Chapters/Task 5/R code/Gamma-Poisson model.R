getwd() # current working directory of R
setwd("/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №5/Chapters/Task 5/R code/")

update_theta <- function(x, sigma, k) {
    theta_values <- array(dim = 3)

    for (i in 1:m) {
        shape <- k + m * mean(x[, i])
        scale <- sigma / (1.0 + sigma * m)

        theta_values[i] <- rgamma(1, shape, scale)
    }

    return(theta_values)
}

update_sigma <- function(x, theta, k) {
    shape <- n * k + 1
    scale <- n * mean(theta) + 1

    sigma_value <- 1 / rgamma(1, shape, scale)

    return(sigma_value)
}

gibbs <- function(data, n_iter, init, prior) {
    # initialize
    sigma_out <- array(data = 0.0, dim = c(n_iter, 1))
    theta_out <- array(data = 0.0, dim = c(n_iter, 3))

    sigma_now <- init$sigma

    # Gibbs sampler
    for (i in 1:n_iter) {
        theta_now <- update_theta(x = data, sigma = sigma_now, k = prior$k)
        sigma_now <- update_sigma(x = data, theta = theta_now, k = prior$k)

        theta_out[i, ] <- theta_now
        sigma_out[i, ] <- sigma_now
    }

    cbind(sigma = sigma_out[, 1], theta1 = theta_out[, 1], theta2 = theta_out[, 2], theta3 = theta_out[, 3])
}

# Load the dataset
data <- read.csv("hw5.csv")
head(data)
summary(data)

x1 <- data$V1
x2 <- data$V2
x3 <- data$V3

hist(x1, freq = FALSE)
hist(x2, freq = FALSE)
hist(x3, freq = FALSE)

x <- cbind(x1, x2, x3)

n <- nrow(x)
m <- ncol(x)

prior <- list()
prior$k <- 3

set.seed(53)

init <- list()
init$sigma <- 1.0

post <- gibbs(data = x, n_iter = 5e3, init = init, prior = prior)
post

library("coda")
plot(as.mcmc(post))

summary(as.mcmc(post))