#' @title RLB Class
#' @description A class to represent a Multiple Linear Regression model.
#' @slot y A numeric vector representing the response variable.
#' @slot X A numeric matrix representing the predictor variables.
#' @slot koef A numeric vector representing the estimated coefficients.
#' @exportClass RLB
setClass("RLB",
         slots = list(
           y = "numeric",
           X = "matrix",
           koef = "numeric"
         ))

#' @title Constructor for RLB Class
#' @description Creates an object of class RLB representing a Multiple Linear Regression model.
#' @param y A numeric vector for the response variable.
#' @param ... One or more numeric vectors or a matrix for the predictor variables.
#' These will be column-binded to form the predictor matrix X.
#' @return An object of class \code{RLB}.
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' x1_data <- rnorm(n)
#' x2_data <- rnorm(n)
#' y_data <- 5 + 2*x1_data - 3*x2_data + rnorm(n)
#' model <- RLB(y_data, x1_data, x2_data)
#' print(model)
#' }
RLB <- function(y, ...){
  nrow <- length(y)
  x <- cbind(...)

  if(!is.numeric(y) || !is.numeric(x)) stop("X dan y harus numerik")
  if(nrow(x) != nrow) stop("Jumlah amatan x dan y harus sama")

  fungsi <- function(params, y, x){
    y.duga <- cbind(1, x) %*% as.matrix(params)
    smsq.resid <- sum((y - y.duga)^2)
    return(smsq.resid)
  }

  hasil.optim <- optim(rep(1, 1 + ncol(x)), fungsi, y = y, x = x)
  new("RLB", y = y, X = x, koef = hasil.optim$par)
}

#' @title Get Response Variable
#' @description Extracts the response variable from an RLB object.
#' @param model An object of class \code{RLB}.
#' @return A numeric vector representing the response variable y.
#' @export
#' @docType methods
#' @rdname respons-methods
respons <- function(model) model@y

#' @title Get Predictor Variables
#' @description Extracts the predictor variables matrix from an RLB object.
#' @param model An object of class \code{RLB}.
#' @return A numeric matrix representing the predictor variables X.
#' @export
#' @docType methods
#' @rdname penjelas-methods
penjelas <- function(model) model@X

#' @title Get Regression Coefficients
#' @description Extracts the estimated regression coefficients from an RLB object.
#' @param model An object of class \code{RLB}.
#' @return A numeric vector representing the estimated coefficients.
#' @export
#' @docType methods
#' @rdname koef.reg-methods
koef.reg <- function(model) model@koef

#' @title Plot Method for RLB Objects
#' @description Generates scatter plots of each predictor variable against the response variable,
#' with the simple linear regression line overlaid.
#' @param x An object of class \code{RLB}.
#' @param y Missing, not used in this method.
#' @param ... Additional arguments passed to the plot function.
#' @importFrom graphics plot abline par
#' @importFrom stats lm
#' @method plot RLB
#' @export
#' @docType methods
#' @rdname plot-RLB-methods
setMethod("plot", signature(x = "RLB", y = "missing"),
          function(x, y, ...) {
            X_plot <- penjelas(x)
            y_plot <- respons(x)

            par(mfrow = c(1, ncol(X_plot)))
            for (i in 1:ncol(X_plot)) {
              plot(X_plot[,i], y_plot, pch=19, col="navy",
                   xlab = paste0("X", i),
                   ylab = "y",
                   main = paste("X", i, "vs y"))
              abline(lm(y_plot ~ X_plot[,i]), col="red")
            }
          })

#' @title Residuals Method for RLB Objects
#' @description Calculates the residuals for an RLB model.
#' @param object An object of class \code{RLB}.
#' @param ... Additional arguments (not used).
#' @return A numeric vector representing the residuals (y - y_hat).
#' @importFrom stats residuals
#' @export
#' @docType methods
#' @rdname residuals-RLB-methods
setMethod("residuals", signature(object = "RLB"),
          function(object, ...){
            y_hat <- cbind(1, object@X) %*% as.matrix(object@koef)
            respons(object) - y_hat
          })

#' @title Summary Method for RLB Objects
#' @description Provides a summary of the RLB model, including coefficients,
#' standard errors, t-statistics, p-values, R-squared, and other metrics.
#' @param object An object of class \code{RLB}.
#' @param ... Additional arguments (not used).
#' @return Prints a summary of the model to the console. Invisibly returns a list
#' containing the coefficient table and other summary statistics (this part is not
#' explicitly implemented in the original code, it only prints).
#' @importFrom stats pt solve
#' @method summary RLB
#' @export
#' @docType methods
#' @rdname summary-RLB-methods
setMethod("summary", signature(object = "RLB"),
          function(object, ...) {
            X <- cbind(1, object@X)
            y <- object@y
            beta <- object@koef

            y_hat <- X %*% as.matrix(beta)
            resid <- y - y_hat
            n <- length(y)
            p <- ncol(X)

            SSE <- sum(resid^2)
            MSE <- SSE / (n - p)
            RMSE <- sqrt(MSE)

            SST <- sum((y - mean(y))^2)
            R2 <- 1 - SSE / SST
            R2_adj <- 1 - ((1 - R2) * (n - 1)) / (n - p)

            XtX_inv <- tryCatch(solve(t(X) %*% X), error = function(e) {
              warning("Could not invert X'X matrix. Standard errors, t-stats, and p-values may be unreliable.")
              matrix(NA, nrow = p, ncol = p)
            })

            se_beta <- sqrt(diag(MSE * XtX_inv))
            t_stat <- beta / se_beta
            p_val <- 2 * (1 - pt(abs(t_stat), df = n - p))

            coef_table <- data.frame(
              Estimate = round(beta, 4),
              Std.Error = round(se_beta, 4),
              t.value = round(t_stat, 4),
              p.value = round(p_val, 4),
              Signif = ifelse(p_val < 0.05, "*", "")
            )
            rownames(coef_table) <- paste0("beta", 0:(p - 1))

            cat("=== Model Regresi Linier Berganda (RLB) ===\n\n")
            cat("Koefisien Estimasi:\n")
            print(round(object@koef, 4))
            print(coef_table)

            cat("\nResidual Standard Error (RMSE):", round(RMSE, 4), "\n")
            cat("R-squared: ", round(R2, 4), "\n")
            cat("Adjusted R-squared: ", round(R2_adj, 4), "\n")
            cat("Mean Squared Error (MSE):", round(MSE, 4), "\n")
            cat("Jumlah Amatan:", n, "\n")
            cat("Jumlah Prediktor:", p - 1, "\n")
            cat("* signifikan di alpha 5%\n")
          })

#' @title RLB_stepwise Class
#' @description A class that inherits from RLB and implements stepwise variable selection
#' for multiple linear regression.
#' @slot y A numeric vector representing the response variable.
#' @slot X A numeric matrix representing the predictor variables selected by stepwise regression.
#' @slot koef A numeric vector representing the estimated coefficients for the selected model.
#' @exportClass RLB_stepwise
setClass("RLB_stepwise", contains = "RLB")

#' @title Stepwise Variable Selection
#' @description Performs stepwise variable selection using AIC.
#' This is an internal helper function.
#' @param X A numeric matrix of potential predictor variables.
#' @param y A numeric vector for the response variable.
#' @return A matrix containing only the selected predictor variables.
#' If no variables are selected, returns a matrix with 0 columns.
#' @importFrom stats lm step coef as.formula
#' @keywords internal
stepwise_selection <- function(X, y) {
  df <- as.data.frame(X)
  df$y <- y

  constant_cols <- sapply(df[, -ncol(df), drop = FALSE], function(col) length(unique(col)) == 1)
  if(any(constant_cols)){
    warning("Constant columns found in predictors and were removed before stepwise selection: ",
            paste(names(df)[constant_cols], collapse=", "))
    df <- df[, !constant_cols]
    if (ncol(df) <= 1) {
      return(matrix(ncol = 0, nrow = nrow(X)))
    }
  }

  if (ncol(df) == 1 && "y" %in% names(df)) {
    return(matrix(ncol = 0, nrow = nrow(X)))
  }

  predictor_names <- setdiff(names(df), "y")
  if (length(predictor_names) == 0) {
    return(matrix(ncol = 0, nrow = nrow(X)))
  }

  formula_str_full <- paste("y ~", paste(predictor_names, collapse = " + "))
  full_model <- lm(as.formula(formula_str_full), data = df)
  null_model <- lm(y ~ 1, data = df)

  step_model <- suppressMessages(step(null_model,
                                      scope = list(lower = null_model, upper = full_model),
                                      direction = "both", trace = 0, k = 2))

  selected_vars_all <- names(coef(step_model))
  selected_vars <- selected_vars_all[selected_vars_all != "(Intercept)"]

  if (length(selected_vars) == 0) {
    return(matrix(ncol = 0, nrow = nrow(X)))
  }

  X_selected <- as.matrix(df[, selected_vars, drop = FALSE])
  colnames(X_selected) <- selected_vars
  return(X_selected)
}

#' @title Constructor for RLB_stepwise Class
#' @description Creates an object of class RLB_stepwise, performing stepwise
#' variable selection before fitting the Multiple Linear Regression model.
#' @param y A numeric vector for the response variable.
#' @param ... One or more numeric vectors or a matrix for the potential predictor variables.
#' These will be column-binded to form the initial predictor matrix.
#' @return An object of class \code{RLB_stepwise}.
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' x1_data <- rnorm(n)
#' x2_data <- rnorm(n)
#' x3_data <- rnorm(n)
#' y_data <- 5 + 2*x1_data + 0*x2_data + 3*x3_data + rnorm(n)
#'
#' model_sw <- RLB_stepwise(y_data, x1_data, x2_data, x3_data)
#' summary(model_sw)
#' plot(model_sw)
#' }
RLB_stepwise <- function(y, ...) {
  x_full <- cbind(...)

  if (is.null(colnames(x_full)) && ncol(x_full) > 0) {
    colnames(x_full) <- paste0("X", 1:ncol(x_full))
  }

  nrow <- length(y)

  if (!is.numeric(y) || !is.numeric(x_full)) stop("X dan y harus numerik")
  if (nrow(x_full) != nrow) stop("Jumlah amatan x dan y harus sama")
  if (ncol(x_full) == 0) stop("Tidak ada variabel prediktor yang diberikan.")


  x_selected <- stepwise_selection(x_full, y)

  fungsi <- function(params, y, x){
    if (ncol(x) == 0) {
      y.duga <- rep(params[1], length(y))
    } else {
      y.duga <- cbind(1, x) %*% as.matrix(params)
    }
    smsq.resid <- sum((y - y.duga)^2)
    return(smsq.resid)
  }

  num_params <- 1 + ncol(x_selected)
  initial_params <- rep(1, num_params)

  hasil.optim <- optim(initial_params, fungsi, y = y, x = x_selected)
  new("RLB_stepwise", y = y, X = x_selected, koef = hasil.optim$par)
}
