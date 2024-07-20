# install.packages("gslnls")
library("ggplot2")
library("patchwork")

dEUc <- function(T, Ub, Lb, Db, Cb) {
    return(
        0.5 + 0.5 * tanh((T - Db) / Cb)
    )
}

dELc <- function(T, Ub, Lb, Db, Cb) {
    return(
        0.5 - 0.5 * tanh((T - Db) / Cb)
    )
}

dEDc <- function(T, Ub, Lb, Db, Cb) {
    return(
        (Ub + Lb) / 2 - ((Ub - Lb) / 2) * (1 / cosh((T - Db) / Cb))
    )
}

dECc <- function(T, Ub, Lb, Db, Cb) {
    return(
        (Ub + Lb) / 2 - ((Ub - Lb) * (T - Db) / (2 * Cb^2)) * (1 / cosh((T - Db) / Cb))
    )
}

Eb <- function(T, Ub, Lb, Db, Cb) {
    return(
        (Ub + Lb) / 2 + ((Ub - Lb) / 2) * tanh((T - Db) / Cb)
    )
}

energy_function <- function(T, U, L, D, C) {
    return(
        (U + L) / 2 + ((U - L) / 2) * tanh((T - D) / C)
    )
}

reverse_energy_function <- function(E, U, L, D, C) {
    return(
        (C * log((E - L) / (U - E))) / 2 + D
    )
}

input_parameters <- data.frame(
    U = 210.0,
    L = 2.0,
    D = -50.0,
    C = 15.0,
    initial_guess <- "random", # "random", "input values"
    points_number <- 12,
    samples_per_point <- 1,
    nls_runs_number <- 20,
    nls_iterations_per_run <- 300
)

Uadm <- 10
Ladm <- 0.25
Dadm <- 0.25
Cadm <- 0.2

t_values <- seq(-90, 0, length.out = input_parameters$points_number)
# t_values <- append(t_values, c(-190, 100)) # attempt to stabilize solution

estimated_coeffs_dataframe <- data.frame()

coef_labels <- c("U", "L", "D", "C")
coef_histplot <- list()
coef_traceplot <- list()

images_dir <- "/home/anton/Code/KPI-FAM-subjects/Weekly diploma reports/Homework №2/LaTeX/Images/"

for (run in 1:input_parameters$nls_runs_number) {
    # Generate sample points
    x_points <- y_points <- array(
        data = 0.0,
        dim = c(input_parameters$points_number, input_parameters$samples_per_point)
    )

    for (i in 1:input_parameters$points_number) {
        for (j in 1:input_parameters$samples_per_point) {
            x_points[i, j] <- t_values[i]
            y_points[i, j] <- energy_function(
                T = t_values[i],
                U = input_parameters$U + rnorm(1, mean = 0, sd = 5),
                L = input_parameters$L + rnorm(1, mean = 0, sd = 5),
                D = input_parameters$D + rnorm(1, mean = 0, sd = 5),
                C = input_parameters$C + rnorm(1, mean = 0, sd = 5)
            )
        }
    }

    data <- data.frame(T = as.vector(x_points), E = as.vector(y_points))

    if (input_parameters$initial_guess == "input values") {
        Ub <- input_parameters$U
        Lb <- input_parameters$L
        Db <- input_parameters$D
        Cb <- input_parameters$C
    } else if (input_parameters$initial_guess == "random") {
        Ub <- runif(1, 50, 150)
        Lb <- runif(1, 5, 15)
        Db <- runif(1, -40, -10)
        Cb <- runif(1, 5, 10)
    }

    iteration <- 1
    diff_error <- 100
    estimated_coeffs_per_run <- data.frame()

    while (iteration <= input_parameters$nls_iterations_per_run && diff_error > 0.1) {
        main_matrix <- array(data = 0.0, dim = c(4, 4))
        free_matrix <- array(data = 0.0, dim = c(4))

        for (i in 1:length(data$T)) {
            main_matrix[1, 1] <- main_matrix[1, 1] + dEUc(data$T[i], Ub, Lb, Db, Cb) * dEUc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[1, 2] <- main_matrix[1, 2] + dELc(data$T[i], Ub, Lb, Db, Cb) * dEUc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[1, 3] <- main_matrix[1, 3] + dEDc(data$T[i], Ub, Lb, Db, Cb) * dEUc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[1, 4] <- main_matrix[1, 4] + dECc(data$T[i], Ub, Lb, Db, Cb) * dEUc(data$T[i], Ub, Lb, Db, Cb)

            main_matrix[2, 1] <- main_matrix[2, 1] + dEUc(data$T[i], Ub, Lb, Db, Cb) * dELc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[2, 2] <- main_matrix[2, 2] + dELc(data$T[i], Ub, Lb, Db, Cb) * dELc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[2, 3] <- main_matrix[2, 3] + dEDc(data$T[i], Ub, Lb, Db, Cb) * dELc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[2, 4] <- main_matrix[2, 4] + dECc(data$T[i], Ub, Lb, Db, Cb) * dELc(data$T[i], Ub, Lb, Db, Cb)

            main_matrix[3, 1] <- main_matrix[3, 1] + dEUc(data$T[i], Ub, Lb, Db, Cb) * dEDc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[3, 2] <- main_matrix[3, 2] + dELc(data$T[i], Ub, Lb, Db, Cb) * dEDc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[3, 3] <- main_matrix[3, 3] + dEDc(data$T[i], Ub, Lb, Db, Cb) * dEDc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[3, 4] <- main_matrix[3, 4] + dECc(data$T[i], Ub, Lb, Db, Cb) * dEDc(data$T[i], Ub, Lb, Db, Cb)

            main_matrix[4, 1] <- main_matrix[4, 1] + dEUc(data$T[i], Ub, Lb, Db, Cb) * dECc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[4, 2] <- main_matrix[4, 2] + dELc(data$T[i], Ub, Lb, Db, Cb) * dECc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[4, 3] <- main_matrix[4, 3] + dEDc(data$T[i], Ub, Lb, Db, Cb) * dECc(data$T[i], Ub, Lb, Db, Cb)
            main_matrix[4, 4] <- main_matrix[4, 4] + dECc(data$T[i], Ub, Lb, Db, Cb) * dECc(data$T[i], Ub, Lb, Db, Cb)

            free_matrix[1] <- free_matrix[1] + (
                data$E[i] * dEUc(data$T[i], Ub, Lb, Db, Cb)
            ) - (
                Eb(data$T[i], Ub, Lb, Db, Cb) * dEUc(data$T[i], Ub, Lb, Db, Cb)
            )

            free_matrix[2] <- free_matrix[2] + (
                data$E[i] * dELc(data$T[i], Ub, Lb, Db, Cb)
            ) - (
                Eb(data$T[i], Ub, Lb, Db, Cb) * dELc(data$T[i], Ub, Lb, Db, Cb)
            )

            free_matrix[3] <- free_matrix[3] + (
                data$E[i] * dEDc(data$T[i], Ub, Lb, Db, Cb)
            ) - (
                Eb(data$T[i], Ub, Lb, Db, Cb) * dEDc(data$T[i], Ub, Lb, Db, Cb)
            )

            free_matrix[4] <- free_matrix[4] + (
                data$E[i] * dECc(data$T[i], Ub, Lb, Db, Cb)
            ) - (
                Eb(data$T[i], Ub, Lb, Db, Cb) * dECc(data$T[i], Ub, Lb, Db, Cb)
            )
        }

        estimated_coeffs <- solve(main_matrix, free_matrix)

        Uc <- estimated_coeffs[1]
        Lc <- estimated_coeffs[2]
        Dc <- estimated_coeffs[3]
        Cc <- estimated_coeffs[4]

        estimated_coeffs_per_run <- rbind(estimated_coeffs_per_run,
            data.frame(
                U = Ub + Uc,
                L = Lb + Lc,
                D = Db + Dc,
                C = Cb + Cc
            )
        )

        kU <- Uc / Uadm
        kL <- Lc / Ladm
        kD <- Dc / Dadm
        kC <- Cc / Cadm
        kmax <- max(kU, kL, kD, kC, 1)

        Ub <- Ub + Uc / kmax
        Lb <- Lb + Lc / kmax
        Db <- Db + Dc / kmax
        Cb <- Cb + Cc / kmax

        diff_error <- max(kU, kL, kD, kC)
        iteration <- iteration + 1
    }

    for (i in 1:4) {
        coef_traceplot[[i]] <- ggplot() +
            geom_point(
                data = data.frame(x = 1:length(estimated_coeffs_per_run[, 1]), y = estimated_coeffs_per_run[, i]),
                mapping = aes(x, y)
            ) +
            geom_hline(
                data = data.frame(y = input_parameters[, i]),
                mapping = aes(yintercept = y, color = "Real value"),
                lwd = 0.8
            ) +
            labs(
                x = "Iteration",
                y = "Estimated value",
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

        ggsave(
            filename = paste(
                images_dir, "NLS manual: ", coef_labels[i], " traceplot (last nls run)", ".png",
                sep = ""
            ),
            plot = coef_traceplot[[i]],
            width = 7,
            height = 5,
            units = "in",
            dpi = 1000,
        )
    }

    estimated_coeffs_dataframe <- rbind(
        estimated_coeffs_dataframe,
        tail(estimated_coeffs_per_run, 1),
        make.row.names = FALSE
    )
}

for (i in 1:4) {
    coef_histplot[[i]] <- ggplot() +
        geom_histogram(
            data = data.frame(x = estimated_coeffs_dataframe[, i][
                estimated_coeffs_dataframe[, i] > -500 & estimated_coeffs_dataframe[, i] < 500
            ]),
            mapping = aes(x, y = after_stat(density)),
            # DENSITY INFO: https://plotnine.org/reference/geom_histogram
            # LEGEND: add fill = "Simulation histogram" inside aes() in order to legend
            fill = "gray",
            color = "#6f6f6f",
            bins = 40
        ) +
        geom_vline(
            data = data.frame(x = input_parameters[, i]),
            mapping = aes(xintercept = x, color = "Real value"),
            lwd = 0.8
        ) +
        labs(
            x = paste(coef_labels[i], "coefficient value"),
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

    ggsave(
        filename = paste(
            images_dir, "NLS manual: ", coef_labels[i], " histplot", ".png",
            sep = ""
        ),
        plot = coef_histplot[[i]],
        width = 7,
        height = 5,
        units = "in",
        dpi = 1000,
    )
}

for (i in 1:4) {
    print(paste(coef_labels[i], "outliers:", sum(
        estimated_coeffs_dataframe[, i] < -500 | estimated_coeffs_dataframe[, i] > 500,
        na.rm = TRUE
    )))
}

true_temperature <- reverse_energy_function(
    E = 110,
    U = input_parameters$U,
    L = input_parameters$L,
    D = input_parameters$D,
    C = input_parameters$C
)

temperature <- array()
for (i in 1:input_parameters$nls_runs_number....10) {
    temperature[i] <- reverse_energy_function(
        E = 110,
        U = estimated_coeffs_dataframe[i, 1],
        L = estimated_coeffs_dataframe[i, 2],
        D = estimated_coeffs_dataframe[i, 3],
        C = estimated_coeffs_dataframe[i, 4]
    )
}

temperature_histplot <- ggplot() +
    geom_histogram(
        data = data.frame(x = temperature),
        mapping = aes(x, y = after_stat(density)),
        # DENSITY INFO: https://plotnine.org/reference/geom_histogram
        # LEGEND: add fill = "Simulation histogram" inside aes() in order to legend
        fill = "gray",
        color = "#6f6f6f",
        bins = 40
    ) +
    geom_vline(
        data = data.frame(x = true_temperature),
        mapping = aes(xintercept = x, color = "Real value"),
        lwd = 0.8
    ) +
    labs(
        x = "Temperature value",
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

ggsave(
    filename = paste(
        images_dir, "NLS manual: temperature histplot.png",
        sep = ""
    ),
    plot = temperature_histplot,
    width = 7,
    height = 5,
    units = "in",
    dpi = 1000,
)

print(paste("Temperature NaNs:", sum(is.na(temperature))))