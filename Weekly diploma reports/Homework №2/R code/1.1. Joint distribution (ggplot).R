# install.packages("dplyr")
library("gslnls")
library("ggplot2")
library("patchwork")

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

plot_hist_and_distribution <- function(t, t0, a, b, mu_c, sigma_c, sample_size, title_name, plot_info) {
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

    # Pull just the histogram info to smoothly draw geom_function() curve
    geom_hist_data <- ggplot_build(
        ggplot() +
            geom_histogram(
                data = data.frame(x = energy[energy_limit_left < energy & energy < energy_limit_right]),
                mapping = aes(x, after_stat(density)),
                bins = 40
            )
    )

    normal_distribution_plot <- ggplot() +
        geom_histogram(
            data = data.frame(x = energy[energy_limit_left < energy & energy < energy_limit_right]),
            mapping = aes(x, after_stat(density)),
            # DENSITY INFO: https://plotnine.org/reference/geom_histogram
            # LEGEND: add fill = "Simulation histogram" inside aes() in order to legend
            fill = "gray",
            color = "#6d6d6d",
            bins = 40,
            position = position_dodge(.7)
        ) +
        geom_function(
            mapping = aes(color = "Analytical distribution"),
            xlim = c(
                min(geom_hist_data$data[[1]]$xmin),
                max(geom_hist_data$data[[1]]$xmax)
            ),
            fun = energy_density,
            n = 500,
            na.rm = TRUE,
            args = list(t, t0, a, b, mu_c, sigma_c),
            # color = "blue",
        ) +
        labs(x = plot_info$xlabel, y = plot_info$ylabel, title = title_name) +
        # scale_fill_manual(name = NULL, values = c("gray")) + # Legend for fill color
        scale_color_manual(name = NULL, values = c("blue")) + # Legend for line color
        theme(
            legend.position = plot_info$legend_position,
            legend.justification = c(1.0, 1.0),
            legend.box.background = element_rect(color = "black", size = 1),
            legend.text = element_text(family = "CMU serif", size = 12),
            axis.text = element_text(family = "CMU serif", size = 10),
            axis.title = element_text(family = "CMU serif", size = 12),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            plot.title = element_text(family = "CMU serif", hjust = 0.5, size = 12)
        )

    return(normal_distribution_plot)
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
mode_index <- 3

if (mode_index != 3) {
    plot_info <- list(
        width = 8,
        height = 5,
        xlabel = "Energy value",
        ylabel = "Density",
        legend_position = c(0.98, 0.98)
    )

    normal_distribution_plot <- plot_hist_and_distribution(
        t = experiment$t_value[mode_index],
        t0 = experiment$t0[mode_index],
        a = experiment$a[mode_index],
        b = experiment$b[mode_index],
        mu_c = experiment$mu_c[mode_index],
        sigma_c = experiment$sigma_c[mode_index],
        sample_size = sample_size,
        title_name = NULL,
        plot_info = plot_info
    )
} else if (mode_index == 3) {
    plot_info <- list(
        width = 14,
        height = 8,
        xlabel = NULL,
        ylabel = NULL,
        legend_position = "none"
    )

    normal_distribution_subplots <- list()
    t_values <- seq(-90, 0, length.out = 12)

    for (i in 1:12) {
        normal_distribution_subplots[[i]] <- plot_hist_and_distribution(
            t = t_values[i],
            t0 = experiment$t0[mode_index],
            a = experiment$a[mode_index],
            b = experiment$b[mode_index],
            mu_c = experiment$mu_c[mode_index],
            sigma_c = experiment$sigma_c[mode_index],
            sample_size = sample_size,
            title_name = paste("T =", round(t_values[i], 4)),
            plot_info = plot_info
        )
    }

    normal_distribution_plot <- wrap_plots(normal_distribution_subplots, ncol = 4, nrow = 3)
}

ggsave(
    filename = paste(
        "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №2/LaTeX/Images/",
        experiment$mode[mode_index], ".png",
        sep = ""
    ),
    plot = normal_distribution_plot,
    width = plot_info$width,
    height = plot_info$height,
    units = "in",
    dpi = 1000,
)