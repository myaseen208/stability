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

    Y   <- deparse(substitute(.y))
    Rep <- deparse(substitute(.rep))
    G   <- deparse(substitute(.gen))
    E   <- deparse(substitute(.env))

    g <- length(levels(.data[[G]]))
    e <- length(levels(.data[[E]]))
    r <- length(levels(.data[[Rep]]))



    g_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G)) %>%
      dplyr::summarize(Mean = mean(!!rlang::sym(Y)))


    ge_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G), !!rlang::sym(E)) %>%
      dplyr::summarize(GE.Mean = mean(!!rlang::sym(Y))) %>%
      tidyr::spread(key = E, value = GE.Mean)

    ge_means1 <- as.matrix(ge_means[, -1])
    rownames(ge_means1) <- c(ge_means[, 1])[[1]]

    gge_effects <-
      sweep(
        x      = ge_means1
        , MARGIN = 2
        , STATS  = colMeans(ge_means1)
      )

    ge_effects <-
      sweep(
        x      = gge_effects
        , MARGIN = 1
        , STATS  = rowMeans(gge_effects)
      )


    SVD <- svd(ge_effects)
    PC <- SVD$u %*% diag(sqrt(SVD$d))
    dimnames(PC) <- list(row.names(ge_effects), paste0("PC", 1:e))

    Lambda  <- SVD$d

    MASD <-
      g_means %>%
      dplyr::mutate(
          MASD  = sqrt(diag(PC[ ,1:(.m)] %*% diag(Lambda[1:(.m)]) %*% t(PC[ ,1:(.m)])))
        , rMean = min_rank(desc(Mean))
        , rMASD = min_rank(MASD)
        , MYSI  = rMean + rMASD
      )

    return(list(
          MASD = MASD
           ))
  }
