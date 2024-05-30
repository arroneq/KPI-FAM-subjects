# install.packages("dplyr")
library("gslnls")
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

energy_function <- function(T, U, L, C, D) {
    return((U + L) / 2 + ((U - L) / 2) * tanh((T - D) / C))
}

input_parameters <- data.frame(
    U = 210.0,
    L = 2.0,
    D = -50.0,
    C = 15.0,
    points_number <- 12,
    samples_per_point <- 1,
    nls_runs_number <- 1,
    nls_iterations_per_run <- 300
)

Ub <- 100.0
Lb <- 2.0
Db <- 4.0
Cb <- -5.0

Uadm <- 10
Ladm <- 0.25
Dadm <- 0.25
Cadm <- 0.2

t_values <- seq(-90, 0, length.out = input_parameters$points_number)
estimated_coeffs_dataframe <- data.frame()

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
                U = input_parameters$U + rnorm(1, mean = 0, sd = 1),
                L = input_parameters$L + rnorm(1, mean = 0, sd = 1),
                C = input_parameters$C + rnorm(1, mean = 0, sd = 1),
                D = input_parameters$D + rnorm(1, mean = 0, sd = 1)
            )
        }
    }

    data <- data.frame(T = as.vector(x_points), E = as.vector(y_points))

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
        kmax <- max(kU, kL, kD, kC)

        Ub <- Ub + Uc / kmax
        Lb <- Lb + Lc / kmax
        Db <- Db + Dc / kmax
        Cb <- Cb + Cc / kmax

        diff_error <- kmax
        iteration <- iteration + 1
    }

    estimated_coeffs_dataframe <- rbind(
        estimated_coeffs_dataframe,
        estimated_coeffs_per_run[input_parameters$nls_iterations_per_run, ],
        make.row.names = FALSE
    )
}

print(estimated_coeffs_dataframe)
