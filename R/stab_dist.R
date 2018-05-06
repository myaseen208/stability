if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
     "Dist"
  , "rDist"
  ))
}

#' @name    stab_dist
#' @aliases stab_dist
#' @title   Additive ANOVA for Genotypes by Environment Interaction (GEI) model
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
#' Muhammad Yaseen (\email{myaseen208@@gmail.com})
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


    stabDist <-
      GMeans %>%
      dplyr::mutate(
          Dist    = sqrt(diag(PC[ ,1:(!!m)] %*% t(PC[ ,1:(!!m)])))
        , rMean   = min_rank(desc(Mean))
        , rDist   = min_rank(Dist)
        , YSIDist = rMean + rDist
      )

    return(list(
        stabDist = stabDist
           ))
  }
