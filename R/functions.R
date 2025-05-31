# ============================
# 1. Kelas Dasar RLB
# ============================
setClass("RLB",
         slots = list(
           y = "numeric",
           X = "matrix",
           koef = "numeric"
         ))

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

respons <- function(model) model@y
penjelas <- function(model) model@X
koef.reg <- function(model) model@koef

# ============================
# 2. Method tambahan untuk RLB
# ============================
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

setMethod("residuals", signature(object = "RLB"),
          function(object){
            y_hat <- cbind(1, object@X) %*% as.matrix(object@koef)
            respons(object) - y_hat
          })

setMethod("summary", signature(object = "RLB"),
          function(object) {
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

            XtX_inv <- solve(t(X) %*% X)
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

# ============================
# 3. Pewarisan: RLB_stepwise
# ============================

# Definisikan class turunan
setClass("RLB_stepwise", contains = "RLB")

stepwise_selection <- function(X, y) {
  df <- as.data.frame(X)
  df$y <- y

  full_model <- lm(y ~ ., data = df)
  null_model <- lm(y ~ 1, data = df)

  step_model <- step(null_model,
                     scope = list(lower = null_model, upper = full_model),
                     direction = "both", trace = 0)

  selected_vars <- names(coef(step_model))[-1]  # hilangkan intercept
  if (length(selected_vars) == 0) {
    return(matrix(ncol = 0, nrow = nrow(X)))  # Tidak ada variabel terpilih
  }

  X_selected <- as.matrix(df[, selected_vars, drop = FALSE])
  colnames(X_selected) <- selected_vars
  return(X_selected)
}

# Konstruktor untuk RLB_stepwise
#' @export
RLB_stepwise <- function(y, ...) {
  x_full <- cbind(...)
  nrow <- length(y)

  if (!is.numeric(y) || !is.numeric(x_full)) stop("X dan y harus numerik")
  if (nrow(x_full) != nrow) stop("Jumlah amatan x dan y harus sama")

  x_selected <- stepwise_selection(x_full, y)

  fungsi <- function(params, y, x){
    y.duga <- cbind(1, x) %*% as.matrix(params)
    smsq.resid <- sum((y - y.duga)^2)
    return(smsq.resid)
  }

  hasil.optim <- optim(rep(1, 1 + ncol(x_selected)), fungsi, y = y, x = x_selected)
  new("RLB_stepwise", y = y, X = x_selected, koef = hasil.optim$par)
}
