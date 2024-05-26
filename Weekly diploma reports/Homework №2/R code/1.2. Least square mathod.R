# install.packages("dplyr")
library("gslnls")
library("ggplot2")
library("patchwork")

energy_function <- function(t, t0, a, b, c) {
    return(a + b * tanh((t - t0) / c))
}

input_parameters <- data.frame(
    a = 106.0,
    b = 104.0,
    t0 = -50.0,
    mu_c = 15.0,
    sigma_c = 5.0,
    points_number <- 12,
    samples_per_point <- NaN,
    nls_runs_number <- 10000
)

ggplots_list <- list()
samples_per_point_list <- c(100)
plot_label <- paste(
    "TRUE ERROR: NLS regression, n =",
    input_parameters$points_number * samples_per_point_list[1]
)

test_stat_info <- data.frame()

t_values <- seq(-90, 0, length.out = input_parameters$points_number)
c_real <- input_parameters$mu_c

for (k in 1:length(samples_per_point_list)) {
    input_parameters$samples_per_point <- samples_per_point_list[k]

    c_estimated <- array(dim = input_parameters$nls_runs_number)
    for (run in 1:input_parameters$nls_runs_number) {
        # Generate sample points
        x_points <- y_points <- array(dim = c(input_parameters$points_number, input_parameters$samples_per_point))
        for (i in 1:input_parameters$points_number) {
            for (j in 1:input_parameters$samples_per_point) {
                x_points[i, j] <- t_values[i]
                y_points[i, j] <- energy_function(
                    t = t_values[i],
                    t0 = input_parameters$t0,
                    a = input_parameters$a,
                    b = input_parameters$b,
                    c = c_real
                    # c = rnorm(1, mean = c_real, sd = input_parameters$sigma_c)
                ) + rnorm(1, mean = 0, sd = input_parameters$sigma_c)
            }
        }

        data <- data.frame(t = as.vector(x_points), y = as.vector(y_points))

        t0_nls <- input_parameters$t0
        a_nls <- input_parameters$a
        b_nls <- input_parameters$b

        # Run nonlinear least-squares model
        model_gslnls <- gsl_nls(
            fn = y ~ energy_function(
                t = t,
                t0 = t0_nls,
                a = a_nls,
                b = b_nls,
                c
            ),
            data = data,
            algorithm = "lmaccel",
            start = c(c = 1)
        )

        # Summary of the fitted model: summary(model_nls)

        # Get the estimated parameter value
        c_estimated[run] <- coef(model_gslnls)
    }

    test_stat_info <- rbind(test_stat_info,
        data.frame(
            sample_size = 12 * input_parameters$samples_per_point,
            c_real = round(c_real, 4),
            mean = round(mean(c_estimated), 4),
            std = round(sd(c_estimated), 4)
        )
    )

    ggplots_list[[k]] <- ggplot() +
        geom_histogram(
            data = data.frame(c = c_estimated),
            mapping = aes(x = c, y = after_stat(density)),
            # DENSITY INFO: https://plotnine.org/reference/geom_histogram
            # LEGEND: add fill = "Simulation histogram" inside aes() in order to legend
            fill = "gray",
            color = "black",
            bins = 40
        ) +
        geom_vline(
            data = data.frame(x = c_real),
            mapping = aes(xintercept = x, color = "Real value"),
            lwd = 0.8
        ) +
        # geom_point(
        #     data = data,
        #     mapping = aes(x = t, y = y)
        # ) +
        # geom_function(
        #     xlim = c(min(data$t), max(data$t)),
        #     fun = energy_function,
        #     n = 500,
        #     na.rm = TRUE,
        #     args = list(t0 = t0, a = a, b = b, c = estimated_params),
        #     colour = "red"
        # ) +
        labs(
            x = "C estimation value",
            y = "Density",
            # title = paste("One-shot sample size n =", 12 * input_parameters$samples_per_point)
        ) +
        # scale_fill_manual(name = NULL, values = c("gray")) + # Legend for fill color
        scale_color_manual(name = NULL, values = c("blue")) + # Legend for line color
        theme(
            legend.position = c(0.98, 0.98),
            legend.justification = c(1.0, 1.0),
            # legend.key.width = unit(3, "cm"),
            legend.box.background = element_rect(color = "black", linewidth = 1),
            legend.text = element_text(family = "CMU serif", size = 12),
            axis.text = element_text(family = "CMU serif", size = 10),
            axis.title = element_text(family = "CMU serif", size = 12),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            plot.title = element_text(
                margin = margin(t = 0, r = 0, b = 10, l = 0),
                family = "CMU serif",
                hjust = 0.5,
                size = 12
            )
        )
}

print(test_stat_info)
final_ggplot <- wrap_plots(ggplots_list, ncol = length(samples_per_point_list), nrow = 1)

ggsave(
    filename = paste(
        "/home/anton/Code/KPI FAM subjects/Weekly diploma reports/Homework №2/LaTeX/Images/",
        plot_label,
        ".png",
        sep = ""
    ),
    plot = final_ggplot,
    width = 7,
    height = 5,
    units = "in",
    dpi = 1000,
)