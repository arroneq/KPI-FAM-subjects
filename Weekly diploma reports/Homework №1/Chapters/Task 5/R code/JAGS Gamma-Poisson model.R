# install.packages("ggmcmc")
library("ggplot2")
library("patchwork")
library("ggmcmc")
library("coda")
library("rjags")

# Current working directory of R
getwd()

# Change working directory of R
setwd("/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/R code/")

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

# Create clear dataset
x <- cbind(x1, x2, x3)

# ------------------------------------------------------------------------
# JAGS GIBBS SAMPLER
# ------------------------------------------------------------------------

# Set main parameters
n <- ncol(x)
m <- nrow(x)

prior <- list()
prior$k <- 3

init <- list()
init$sigma <- 1.0

jags_model_string <- "model {
    for (i in 1:n * m) {
        x[i] ~ dpois(theta[group[i]])
    }

    for (j in 1:n) {
        theta[j] ~ dgamma(k, invsigma)
    }

    invsigma ~ dexp(1.0)
}"

set.seed(113)

# Convert three flows of data (x1,x2,x3) into one list and label each datapoint accordingly (group 1.0, 2.0 or 3.0)
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
    k = prior$k
)

model <- jags.model(
    textConnection(jags_model_string),
    data = data_jags,
    n.chains = 3
)

update(model, 100) # burn-in period

params <- c("theta", "invsigma")

model_simulation <- coda.samples(
    model = model,
    variable.names = params,
    n.iter = 10e3
)

gelman.diag(model_simulation)
# gelman.plot(model_simulation)
autocorr.diag(model_simulation)
# autocorr.plot(model_simulation)
effectiveSize(model_simulation[1])

# Compute DIC
dic <- dic.samples(model, n.iter = 1e3)

summary(as.mcmc(model_simulation[1]))

# ------------------------------------------------------------------------
# RESULTS OVERVIEW
# ------------------------------------------------------------------------

# LaTeX project has another folder, so set full path and download image [R-base plot]
png(
    filename = "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/Images/JAGS Gibbs traceplot.png",
    width = 8,
    height = 6,
    units = "in",
    bg = "white",
    res = 700,
)

plot(model_simulation[1][, c(1, 2)])
dev.off()

# In order to display greek letters in a ggplot() correctly
for (chain in 1:3) {
    colnames(model_simulation[[chain]]) <- c("sigma^{-1}", "theta[1]", "theta[2]", "theta[3]")
}

ggmcmc_gelman_diag_plot <- ggs_grb(ggs(model_simulation[, c(1, 2)]), greek = TRUE) +
    theme(
        axis.text = element_text(family = "CMU serif", size = 6),
        axis.title = element_text(family = "CMU serif", size = 8),
        # axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(family = "CMU serif", hjust = 0.5, size = 8)
    )

ggsave(
    filename = "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/Images/JAGS Gibbs gelman diagnostic.png",
    plot = ggmcmc_gelman_diag_plot,
    width = 6,
    height = 3,
    units = "in",
    dpi = 1000,
)

ggmcmc_plots_list <- list()
ggmcmc_plots_list[[1]] <- ggs_traceplot(ggs(model_simulation[, c(1, 2)]), greek = TRUE)
ggmcmc_plots_list[[2]] <- ggs_density(ggs(model_simulation[, c(1, 2)]), greek = TRUE)

ggmcmc_plot <- wrap_plots(ggmcmc_plots_list, ncol = 2, nrow = 1)

ggsave(
    filename = "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/Images/JAGS Gibbs traceplot.png",
    plot = ggmcmc_plot,
    width = 10,
    height = 5,
    units = "in",
    dpi = 1000,
)

ggmcmc_autocorr_plot <- ggs_autocorrelation(ggs(model_simulation[1][, c(1, 2)]), greek = TRUE) +
    theme(
        axis.text = element_text(family = "CMU serif", size = 6),
        axis.title = element_text(family = "CMU serif", size = 8),
        # axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        # axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(family = "CMU serif", hjust = 0.5, size = 12)
    )

ggsave(
    filename = "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №1/Chapters/Task 5/Images/JAGS Gibbs autocorrelation.png",
    plot = ggmcmc_autocorr_plot,
    width = 6,
    height = 3,
    units = "in",
    dpi = 1000,
)

# Calculate Bayesian estimations

theta_mean_estimation <- mean(as.mcmc(model_simulation[1][, 2]))
print(paste("Bayesian mean estimation:", theta_mean_estimation))

x1_mean <- mean(x1)
x_mean <- (sum(x1) + sum(x2) + sum(x3)) / (n * m)

theta_star_estimation <- (prior$k / m + x1_mean) * (x_mean / (prior$k / m + x_mean))
print(paste("Theta star estimation:", theta_star_estimation))