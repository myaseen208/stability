if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
    "GRank"
  ))
}

#' @name    stab_fox
#' @aliases stab_fox
#' @title   Stability Fox Function
#' @description Performs a stability analysis based on the criteria of
#'              Fox et al. (1990), using the statistical "TOP third" only.
#'              In Fox function, a stratified ranking of the genotypes at
#'              each environment separately is done. The proportion of locations
#'              at which the genotype occurred in the top third are expressed in
#'              TOP output.
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
#'          \item  Fox, P.N. and Skovmand, B. and  Thompson, B.K. and  Braun, H.J. and Cormier, R. (1990). Yield and adaptation of hexaploid spring triticale. \emph{Euphytica},  \strong{47}, 57-64.
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
#' YieldFox <-
#'      stab_fox(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'       )
#' YieldFox
#'
stab_fox <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("stab_fox")
}

#' @export
#' @rdname stab_fox

stab_fox.default <-
  function(.data, .y, .rep, .gen, .env){

    Y   <- enquo(.y)
    Rep <- enquo(.rep)
    G   <- enquo(.gen)
    E   <- enquo(.env)

    g <- length(levels(.data$G))
    e <- length(levels(.data$E))
    r <- length(levels(.data$Rep))

    foxOut <-
        .data %>%
          dplyr::group_by(!! E) %>%
          dplyr::mutate(GRank = min_rank(desc(!! Y))) %>%
          dplyr::group_by(!! G) %>%
          dplyr::summarise(
               Mean = mean(!! Y)
            ,  TOP  = sum(GRank <= 3)
            )

    return(list(
       Fox = foxOut
    ))

  }
