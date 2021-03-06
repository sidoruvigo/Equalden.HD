#' Package ‘Equalden.HD’
#'
#' Documentation for package ‘Equalden.HD’ version 1.0
#'
#' @description
#' This package implements three different methods to test the null hypothesis that a large number k of
#' samples have a common density. The sample size can be as small as 2. These methods are particularly
#' well suited to the low sample size, high dimensional setting (n << k). The first method, proposed by
#' Zhan and Hart (2012), was developed to test the null hypothesis when the samples are independent of
#' each other. The other tests, proposed by Cousido-Rocha et al. (2018), are adaptations of the test in Zhan
#' and Hart (2012) for the setting in which the samples are weakly dependent. The standarized version of
#' each test statistic and its p-value are computed among other things.
#'
#'
#' @details
#' \itemize{
#' \item{Package: Equalden.HD}
#' \item{Version: 1.0}
#' \item{Maintainer: Marta Cousido Rocha \email{martacousido@@uvigo.es}}
#' \item{License: GPL-2}
#' }
#'
#' @return
#' \itemize{
#' \item{Equalden.test.HD: Performs the k-sample test proposed in Zhan and Hart (2012) for the setting of low sample size, large
#' dimension and independent samples, and its adaptions to dependent samples proposed in Cousido-Rocha
#' et. al (2018).}
#' }
#'
#' @author
#' \itemize{
#' \item{Cousido Rocha, Marta.}
#' \item{Soage González, José Carlos.}
#' \item{de Uña-Álvarez, Jacobo.}
#' \item{D. Hart, Jeffrey.}
#' }
#'
#' @references
#' \itemize{
#' \item{Cousido-Rocha, M., de Uña-Álvarez, J., and Hart, J.(2018). Testing equality of a large number of densities under mixing conditions. Preprint.}
#' \item{Zhan, D., Hart, J. (2012). Testing equality of a large number of densities. Biometrika, 99, 1-17.}
#' }
#'
#'
"_PACKAGE"
#> [1] "_PACKAGE"
