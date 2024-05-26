# Current working directory of R
getwd()

# Change working directory of R
setwd("/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №2/")

inverse_function <- function(y, t, t0, a, b) {
    return(
        2 * (t - t0) / (log(b - a + y) - log(b + a - y))
    )
}

jacobian <- function(y, t, t0, a, b) {
    return(
        -2 * (t - t0) * (log(b - a + y) - log(b + a - y))^(-2) * (2 * b) / ((b - a + y) * (b + a - y))
    )
}

energy_density <- function(y, t, t0, a, b, mu_c, sigma_c) {
    return(
        dnorm(inverse_function(y, t, t0, a, b), mu_c, sigma_c) * abs(jacobian(y, t, t0, a, b))
    )
}

plot_hist_and_distribution <- function(t, t0, a, b, mu_c, sigma_c, sample_size, title_name) {
    # MCMC generation
    c <- rnorm(sample_size, mu_c, sigma_c)
    energy <-  a + b * tanh((t - t0) / c)

    # Get hist() output in order to exclude short bins
    histogram_output <- hist(
        energy,
        breaks = "Freedman-Diaconis",
        plot = FALSE,
    )

    # Set minimum count in a bin (below is considered as a short bin)
    min_count <- 30

    # Define the smallest count bin in the "left" and "right" sides of a bell-shaped distribution
    energy_leftside_counts <- histogram_output$counts[
        1:which.max(histogram_output$counts)
    ]
    energy_rightside_counts <- histogram_output$counts[
        which.max(histogram_output$counts):length(histogram_output$counts)
    ]

    # Define the the "left" and "right" edges of a bell-shaped distribution
    energy_limit_left <- histogram_output$breaks[
        which(energy_leftside_counts == min(energy_leftside_counts[energy_leftside_counts > min_count]))
    ]
    energy_limit_right <- histogram_output$breaks[
        length(energy_leftside_counts) + which(
            energy_rightside_counts == min(energy_rightside_counts[energy_rightside_counts > min_count])
        )
    ]

    # which() could return two or more values, but we need only one
    energy_limit_left <- head(energy_limit_left, 1)
    energy_limit_right <- tail(energy_limit_right, 1)

    # Simulation histogram
    hist_object <- hist(
        energy[energy_limit_left < energy & energy < energy_limit_right],
        freq = FALSE,
        breaks = "Freedman-Diaconis",
        main = title_name,
        xlab = "Energy value",
        ylab = "Density",
    )

    # Analytical distribution
    curve(
        expr = energy_density(x, t, t0, a, b, mu_c, sigma_c),
        from = min(hist_object$breaks),
        to = max(hist_object$breaks),
        n = 700,
        add = TRUE,
        col = "blue"
    )
}

experiment <- data.frame(
    mode = c(
        "Test №1: naive overview",
        "Test №2: naive phisical contex",
        "Test №3: phisical contex for each T-value"
    ),
    a       = c(1.0, 46.4, 106.0),
    b       = c(1.0, 51.5, 104.0),
    t0      = c(1.0, -64.5, -50.0),
    t_value = c(2.0, -46.0, NaN),
    mu_c    = c(0.0, 18.9, 15.0),
    sigma_c = c(5.0, 5.0, 5.0)
)

sample_size <- 10000
mode_index <- 1

# LaTeX project has another folder, so set full path and download image
png(
    filename = paste(
        "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №2/LaTeX/Images/",
        experiment$mode[mode_index], ".png",
        sep = ""
    ),
    width = 14,
    height = 8,
    units = "in",
    bg = "white",
    res = 700,
)

if (mode_index != 3) {
    plot_hist_and_distribution(
        t = experiment$t_value[mode_index],
        t0 = experiment$t0[mode_index],
        a = experiment$a[mode_index],
        b = experiment$b[mode_index],
        mu_c = experiment$mu_c[mode_index],
        sigma_c = experiment$sigma_c[mode_index],
        sample_size = sample_size,
        title_name = NULL
    )

    legend(
        "topright",
        legend = c("Simulation histogram", "Analytical distribution"),
        pch = c(15, NA),
        col = c("lightgray", "blue"),
        lty = c(NA, 1),
        lwd = c(NA, 2)
    )
} else if (mode_index == 3) {
    par(mfrow = c(3, 4), mar = c(5, 4, 4, 2) + 0.1)

    for (t in seq(-90, 0, length.out = 12)) {
        plot_hist_and_distribution(
            t = t,
            t0 = experiment$t0[mode_index],
            a = experiment$a[mode_index],
            b = experiment$b[mode_index],
            mu_c = experiment$mu_c[mode_index],
            sigma_c = experiment$sigma_c[mode_index],
            sample_size = sample_size,
            title_name = paste("T =", round(t, 4))
        )
    }
}

dev.off()