# install.packages("ggmcmc")

library("ggplot2")
# https://intro2r.com/tips.html#moving-the-legend
# https://r-charts.com/ggplot2/legend/
# https://afit-r.github.io/histograms
# https://ggplot2.tidyverse.org/articles/extending-ggplot2.html

library("patchwork")
library("ggmcmc")
# https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html
# ttps://yutannihilation.github.io/allYourFigureAreBelongToUs/ggmcmc/

library("coda")
library("rjags")
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119942412.app1

# Current working directory of R
getwd()

# Change working directory of R
setwd("/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/R code/")

# ------------------------------------------------------------------------
# FUNCTION DECLARATION
# ------------------------------------------------------------------------

update_theta <- function(x, sigma, k, n, m) {
    theta_values <- array(dim = n)

    for (i in 1:n) {
        shape_theta <- k + m * mean(x[, i])
        scale_theta <- sigma / (1.0 + sigma * m)

        theta_values[i] <- rgamma(1, shape = shape_theta, scale = scale_theta)
    }

    return(theta_values)
}

update_sigma <- function(x, theta, k, n, m) {
    shape_sigma <- n * k + 1
    scale_sigma <- n * mean(theta) + 1

    sigma_value <- 1 / rgamma(1, shape = shape_sigma, scale = 1.0 / scale_sigma)

    return(sigma_value)
}

gibbs <- function(data, n_iter, init, prior, n, m) {
    # Initialize arrays to store variables
    sigma_out <- array(dim = c(n_iter, 1))
    theta_out <- array(dim = c(n_iter, 3))

    sigma_now <- init$sigma

    # Gibbs sampler
    for (i in 1:n_iter) {
        theta_now <- update_theta(x = data, sigma = sigma_now, k = prior$k, n = n, m = m)
        sigma_now <- update_sigma(x = data, theta = theta_now, k = prior$k, n = n, m = m)

        theta_out[i, ] <- theta_now
        sigma_out[i, ] <- sigma_now
    }

    cbind(sigma = sigma_out[, 1], theta1 = theta_out[, 1], theta2 = theta_out[, 2], theta3 = theta_out[, 3])
}

# ------------------------------------------------------------------------
# DATA PREPROCESSING
# ------------------------------------------------------------------------

# Load the dataset
data <- read.csv("hw5.csv")
head(data)
summary(data)

# Extract the data
x1 <- data$V1
x2 <- data$V2
x3 <- data$V3

# Recreate the dataset
x <- cbind(x1, x2, x3)

# ------------------------------------------------------------------------
# GIBBS SAMPLER
# ------------------------------------------------------------------------

# Set main parameters
n <- ncol(x)
m <- nrow(x)

prior <- list()
prior$k <- 3

set.seed(53)

init <- list()
init$sigma <- 1.0

posterior <- gibbs(data = x, n_iter = 10e3 + 100, init = init, prior = prior, n = n, m = m)

# Exclude first 100 iterations (burn-in period)
posterior <- tail(posterior, -100)

# Show just few rows to check if a dataset (for sigma and theta1) is built properly
head(posterior[, c(1, 2)])
tail(posterior[, c(1, 2)])

# ------------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS: R-BASE TOOLS
# ------------------------------------------------------------------------

summary(as.mcmc(posterior))

autocorr.diag(as.mcmc(posterior))
# autocorr.plot(as.mcmc(posterior))
effectiveSize(as.mcmc(posterior))

# LaTeX project has another folder, so set full path and download image [R-base plot]
png(
    filename = "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/Images/Gibbs traceplot.png",
    width = 8,
    height = 6,
    units = "in",
    bg = "white",
    res = 500,
)

plot(as.mcmc(posterior[, c(1, 2)]))
dev.off()

# ------------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS: GGMCMC TOOLS
# ------------------------------------------------------------------------

# In order to display greek letters in a ggplot() correctly
colnames(posterior) <- c("sigma", "theta[1]", "theta[2]", "theta[3]")

ggmcmc_plots_list <- list()
ggmcmc_plots_list[[1]] <- ggs_traceplot(ggs(as.mcmc(posterior[, c(1, 2)])), greek = TRUE)
ggmcmc_plots_list[[2]] <- ggs_density(ggs(as.mcmc(posterior[, c(1, 2)])), greek = TRUE)

ggmcmc_plot <- wrap_plots(ggmcmc_plots_list, ncol = 2, nrow = 1)

ggsave(
    filename = "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/Images/Gibbs traceplot.png",
    plot = ggmcmc_plot,
    width = 8,
    height = 5,
    units = "in",
    dpi = 1000,
)

ggmcmc_autocorr_plot <- ggs_autocorrelation(ggs(as.mcmc(posterior[, c(1, 2)])), greek = TRUE) +
    theme(
        axis.text = element_text(family = "CMU serif", size = 6),
        axis.title = element_text(family = "CMU serif", size = 8),
        # axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(family = "CMU serif", hjust = 0.5, size = 12)
    )

ggsave(
    filename = "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/Images/Gibbs autocorrelation.png",
    plot = ggmcmc_autocorr_plot,
    width = 6,
    height = 3,
    units = "in",
    dpi = 1000,
)

# ------------------------------------------------------------------------
# BAYESIAN ESTIMATIONS
# ------------------------------------------------------------------------

theta_mean_estimation <- mean(posterior[, 2])
print(paste("Bayesian empirical mean estimation:", theta_mean_estimation))

x1_mean <- mean(x1)
x_mean <- (sum(x1) + sum(x2) + sum(x3)) / (n * m)

theta_star_estimation <- (prior$k / m + x1_mean) * (x_mean / (prior$k / m + x_mean))
print(paste("Point estimation:", theta_star_estimation))