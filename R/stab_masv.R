if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
     "rMASD"
  ))
}

#' @name    stab_masv
#' @aliases stab_masv
#' @title   Modified Additive Main Effects and Multiplicative Interacion Stability Value
#' @description Additive ANOVA for Genotypes by Environment Interaction (GEI) model
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .rep   Replication Factor
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#' @param .m     No of PCs retained
#'
#' @return Additive ANOVA
#'
#' @author
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Kent M. Edkridge (\email{keskridge1@@unl.edu})
#'          }
#'
#'
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' data(ge_data)
#' YieldMASV <-
#'      stab_masv(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'           , .m    = 2
#'       )
#' YieldMASV
#'
stab_masv <- function(.data, .y, .rep, .gen, .env, .m = 2) {
  UseMethod("stab_masv")
}

#' @export
#' @rdname stab_masv

stab_masv.default <-
  function(.data, .y, .rep, .gen, .env, .m = 2){

    Y   <- enquo(.y)
    Rep <- enquo(.rep)
    G   <- enquo(.gen)
    E   <- enquo(.env)
    m   <- enquo(.m)

    g <- length(levels(.data$G))
    e <- length(levels(.data$E))
    r <- length(levels(.data$Rep))

    GMeans <-
      ge_means(
         .data  = .data
        , .y    = !! Y
        , .gen  = !! G
        , .env  = !! E
      )$g_means

    names(GMeans) <- c("Gen", "Mean")


    GE.Effs <-
      ge_effects(
        .data  = .data
        , .y    = !! Y
        , .gen  = !! G
        , .env  = !! E
      )$ge_effects

    SVD <- svd(GE.Effs)
    PC <- SVD$u %*% diag(sqrt(SVD$d))
    dimnames(PC) <- list(row.names(GE.Effs), paste0("PC", 1:e))

    Lambda  <- SVD$d

    MASD <-
      GMeans %>%
      dplyr::mutate(
          MASD  = sqrt(diag(PC[ ,1:(!!m)] %*% diag(Lambda[1:(!!m)]) %*% t(PC[ ,1:(!!m)])))
        , rMean = min_rank(desc(Mean))
        , rMASD = min_rank(MASD)
        , MYSI  = rMean + rMASD
      )

    return(list(
          MASD = MASD
           ))
  }
