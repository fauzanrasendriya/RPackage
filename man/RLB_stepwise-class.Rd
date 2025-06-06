% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\docType{class}
\name{RLB_stepwise-class}
\alias{RLB_stepwise-class}
\title{Kelas Turunan RLB dengan Seleksi Fitur Stepwise}
\description{
Mendefinisikan kelas S4 bernama "RLB_stepwise" yang mewarisi dari "RLB".
Kelas ini digunakan untuk model regresi linier berganda dimana variabel penjelas
dipilih menggunakan metode seleksi stepwise (berdasarkan AIC).
}
\section{Slots}{

\describe{
\item{\code{y}}{\code{numeric} Vektor respons.}

\item{\code{X}}{\code{matrix} Matriks variabel penjelas yang telah diseleksi.}

\item{\code{koef}}{\code{numeric} Vektor koefisien regresi untuk variabel terpilih.}
}}

\seealso{
\code{\link{RLB-class}}
}
