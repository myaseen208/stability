if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
     "rMean"
  , "rShukaVar"
  ))
}

#' @name    stab_kang
#' @aliases stab_kang
#' @title   Stability Kang Function
#' @description Performs a stability analysis based on the Kang (1988)
#'              criteria. Kang nonparametric stability (ranksum) uses
#'              both "trait single value" and stability variance (Shukla, 1972),
#'              and the genotype with the lowest ranksum is commonly the most favorable one.
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .rep   Replication Factor
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#'
#' @author
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Kent M. Edkridge (\email{keskridge1@@unl.edu})
#'          }
#'
#'
#'
#' @references
#' \enumerate{
#'          \item  Kang, M.S.  (1988). A rank-sum method for selecting high-yielding, stable corn genotypes. \emph{Cereal Research Communications},  \strong{16}, 1-2.
#'          \item  Shukla, G.K.  (1972). Some aspects of partitioning genotype environmental components of variability. \emph{Heredity},  \strong{29}, 237-245.
#'
#'          }
#'
#' @import dplyr
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' data(ge_data)
#' YieldKang <-
#'      stab_kang(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'       )
#' YieldKang
#'
stab_kang <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("stab_kang")
}

#' @export
#' @rdname stab_kang

stab_kang.default <-
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

    g <- nrow(GE.Effs)
    e <- ncol(GE.Effs)

    Wi <- as.matrix(diag(GE.Effs %*% t(GE.Effs)))
    colnames(Wi) <- "Wi"

    ShuklaVar <- (g*(g-1)*Wi - sum(Wi))/((e-1)*(g-1)*(g-2))
    colnames(ShuklaVar) <- "ShuklaVar"

    Kang <-
      tibble::as_tibble(data.frame(GMeans, ShuklaVar)) %>%
      dplyr::mutate(
          rMean     = min_rank(desc(Mean))
        , rShukaVar = min_rank(ShuklaVar)
        , rStab     = rMean + rShukaVar
      )

    return(list(
           Kang = Kang
           ))
  }
