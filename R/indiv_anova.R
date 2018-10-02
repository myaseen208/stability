#' @name    indiv_anova
#' @aliases indiv_anova
#' @title  Individual ANOVA for Each Environment
#' @description Individual ANOVA for Each Environment
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
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Kent M. Edkridge (\email{keskridge1@@unl.edu})
#'          \item Ghulam Murtaza (\email{gmurtaza208@@gmail.com})
#'          }
#'
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#' @import rlang
#' @import tibble
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.indiv_anova <-
#'          indiv_anova(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'           )
#' Yield.indiv_anova
#'
#'
if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c(
        "."
    )
  )
}

indiv_anova <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("indiv_anova")
}

#' @export
#' @rdname indiv_anova

indiv_anova.default <-
  function(.data, .y, .rep, .gen, .env){


    Y   <- deparse(substitute(.y))
    Rep <- deparse(substitute(.rep))
    G   <- deparse(substitute(.gen))
    E   <- deparse(substitute(.env))

    g <- length(levels(.data[[G]]))
    e <- length(levels(.data[[E]]))
    r <- length(levels(.data[[Rep]]))


    ind_aov <-
      .data %>%
      dplyr::group_by(!!rlang::sym(E)) %>%
      dplyr::do(m1 = summary(aov(!!rlang::sym(Y) ~ !!rlang::sym(G), data = .)))

    attr(ind_aov$m1, "names") <-
        paste0(
          "Analysis of Variance for "
          , levels(.data[[E]])
        )

    return(ind_aov = ind_aov$m1)
  }

