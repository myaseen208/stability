if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
    "GRank"
  ))
}

#' @name    stab_fox
#' @aliases stab_fox
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
