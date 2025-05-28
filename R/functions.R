#' @importFrom methods setClass setMethod new signature
#' @importFrom stats optim lm pt residuals
# Penjelasan: 'summary' generic akan ditemukan karena 'stats' ada di Imports (DESCRIPTION)
#' @importFrom graphics par plot abline
NULL

#' Kelas RLB untuk Regresi Linier Berganda
#' @description Mendefinisikan slot untuk objek RLB.
#' @exportClass RLB
setClass("RLB",
         slots = list(
           y = "numeric",
           X = "matrix",
           koef = "numeric"
         ))

#' Konstruktor untuk Objek RLB
#' @description Fungsi untuk membuat dan mengestimasi model RLB.
#' @param y Variabel dependen (numerik).
#' @param ... Variabel independen (numerik), akan di-cbind.
#' @return Objek kelas RLB.
#' @export
RLB <- function(y, ...){
  nrow <- length(y)
  x <- cbind(...)

  if(!is.numeric(y) & !is.numeric(x)) stop("X dan y harus numerik")
  if(nrow(x) != nrow) stop("Jumlah amatan x dan y harus sama")

  fungsi <- function(params, y, x){
    y.duga <- cbind(1, x) %*% as.matrix(params)
    smsq.resid <- sum((y - y.duga)^2)
    return(smsq.resid)}

  hasil.optim <- optim(rep(1, 1 + ncol(x)), fungsi, y=y, x=cbind(...))
  new("RLB", y = y, X = x, koef = hasil.optim$par)
}

#' Aksesor untuk Variabel Respons (y)
#' @description Mengambil slot y dari objek RLB.
#' @param model Objek kelas RLB.
#' @return Mengembalikan slot y.
#' @export
respons <- function(model) model@y

#' Aksesor untuk Matriks Penjelas (X)
#' @description Mengambil slot X dari objek RLB.
#' @param model Objek kelas RLB.
#' @return Mengembalikan slot X.
#' @export
penjelas <- function(model) model@X

#' Aksesor untuk Koefisien Regresi
#' @description Mengambil slot koef dari objek RLB.
#' @param model Objek kelas RLB.
#' @return Mengembalikan slot koef.
#' @export
koef.reg <- function(model) model@koef

#' Metode Plot untuk Objek RLB
#' @description Membuat plot sebar untuk setiap variabel X terhadap y.
#' @param x Objek kelas RLB.
#' @param y Parameter `y` diabaikan untuk method ini.
#' @param ... Argumen tambahan untuk plot.
#' @export
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
              abline(lm(y_plot ~ X_plot[,i]), col="red")}
          })

#' Metode Residuals untuk Objek RLB
#' @description Menghitung residual dari model RLB.
#' @param object Objek kelas RLB.
#' @return Residual dari model.
#' @export
setMethod("residuals", signature(object = "RLB"),
          function(object){
            y_hat <- cbind(1, object@X) %*% as.matrix(object@koef)
            respons(object) - y_hat
          })

#' Metode Summary untuk Objek RLB
#' @description Menampilkan ringkasan statistik model RLB.
#' @param object Objek kelas RLB.
#' @return Mencetak summary dan tidak mengembalikan nilai eksplisit.
#' @export
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
