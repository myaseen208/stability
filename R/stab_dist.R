if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
     "Dist"
  , "rDist"
  ))
}

#' @name    stab_dist
#' @aliases stab_dist
#' @title   Stability Distance in AMMI
#' @description Stability Distance of Genotypes in Additive ANOVA for Genotypes by Environment Interaction (GEI) model
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .rep   Replication Factor
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#' @param .m     No of PCs retained
#'
#' @return Stability Distance
#'
#' @author
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Kent M. Edkridge (\email{keskridge1@@unl.edu})
#'          }
#'
#'
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' data(ge_data)
#' YieldDist <-
#'      stab_dist(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'           , .m    = 2
#'       )
#' YieldDist
#'
stab_dist <- function(.data, .y, .rep, .gen, .env, .m = 2) {
  UseMethod("stab_dist")
}

#' @export
#' @rdname stab_dist

stab_dist.default <-
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


    stabDist <-
      g_means %>%
      dplyr::mutate(
          Dist    = sqrt(diag(PC[ ,1:(.m)] %*% t(PC[ ,1:(.m)])))
        , rMean   = min_rank(desc(Mean))
        , rDist   = min_rank(Dist)
        , YSIDist = rMean + rDist
      )

    return(list(
      stabDist = stabDist
           ))
  }
