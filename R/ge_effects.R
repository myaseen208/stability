#' @name    ge_effects
#' @aliases ge_effects
#' @title Genotype by Environment Interaction Effects
#' @description Calcuates Genotype by Environment Interaction Effects
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Effects
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
#' @import rlang
#' @import tidyr
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.Effects <-
#'               ge_effects(
#'                 .data  = ge_data
#'                , .y    = Yield
#'                , .gen  = Gen
#'                , .env  = Env
#'                )
#' names(Yield.Effects)
#'
#' Yield.Effects$ge_means
#' Yield.Effects$ge_effects
#' Yield.Effects$gge_effects
#'
#'

ge_effects <- function(.data, .y, .gen, .env) {
  UseMethod("ge_effects")
}


#' @export
#' @rdname ge_effects

ge_effects.default <-
  function(.data, .y, .gen, .env){

    Y   <- enquo(.y)
    G   <- enquo(.gen)
    E   <- enquo(.env)


    ge_means <-
      .data %>%
      dplyr::group_by(!! G, !! E) %>%
      dplyr::summarize(GE.Mean = mean(!! Y)) %>%
      tidyr::spread(key = !! E, value = GE.Mean)

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

    return(list(
          "ge_means"    = ge_means
        , "ge_effects"  = ge_effects
        , "gge_effects" = gge_effects
        ))
  }
