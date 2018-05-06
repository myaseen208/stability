if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
     "rASD"
  ))
}

#' @name    stab_asv
#' @aliases stab_asv
#' @title   Additive ANOVA for Genotypes by Environment Interaction (GEI) model
#' @description Additive ANOVA for Genotypes by Environment Interaction (GEI) model
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .rep   Replication Factor
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
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
#' YieldASV <-
#'      stab_asv(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'       )
#' YieldASV
#'
stab_asv <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("stab_asv")
}

#' @export
#' @rdname stab_asv

stab_asv.default <-
  function(.data, .y, .rep, .gen, .env){

    Y   <- enquo(.y)
    Rep <- enquo(.rep)
    G   <- enquo(.gen)
    E   <- enquo(.env)

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

    ASD <-
      GMeans %>%
      dplyr::mutate(
          ASD   = sqrt((Lambda[1]/Lambda[2])*PC[, 1, drop = FALSE]^2 + PC[, 2, drop = FALSE]^2)
        , rMean = min_rank(desc(Mean))
        , rASD  = min_rank(ASD)
        , YSI   = rMean + rASD
      )

    return(list(
          ASD = ASD
           ))
  }
