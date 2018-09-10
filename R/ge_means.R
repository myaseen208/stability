#' @name    ge_means
#' @aliases ge_means
#' @title Genotype by Environment Interaction Means and Ranks
#' @description Calcuates Genotype by Environment Interaction Means along with their Ranks
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Means and Ranks
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
#' @import rlang
#' @importFrom dplyr group_by summarise
#' @importFrom reshape2 acast
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#'
#' Yield.ge_means <-
#'           ge_means(
#'                 .data  = ge_data
#'                , .y    = Yield
#'                , .gen  = Gen
#'                , .env  = Env
#'                )
#'
#' Yield.ge_means$ge_means
#' Yield.ge_means$ge_ranks
#' Yield.ge_means$g_means
#' Yield.ge_means$e_means
#'
#'

ge_means <- function(.data, .y, .gen, .env) {
  UseMethod("ge_means")
}

#' @export
#' @rdname ge_means

ge_means.default <-
  function(.data, .y, .gen, .env){

    Y   <- enquo(.y)
    G   <- enquo(.gen)
    E   <- enquo(.env)

    ge_means <-
      .data %>%
      dplyr::group_by(!! G, !! E) %>%
      dplyr::summarize(Mean = mean(!! Y)) %>%
      tidyr::spread(key = !! E, value = Mean)

    g_means <-
      .data %>%
      dplyr::group_by(!! G) %>%
      dplyr::summarize(Mean = mean(!! Y))

    e_means <-
      .data %>%
        dplyr::group_by(!! E) %>%
        dplyr::summarize(Mean = mean(!! Y))


    ge_means1 <- as.matrix(ge_means[, -1])
    rownames(ge_means1) <- c(ge_means[, 1])$G

    Eg_means  <- t(ge_means1)
    ge_ranks  <- matrix(data = NA, nrow = nrow(Eg_means), ncol = ncol(Eg_means))
    for(i in 1:nrow(Eg_means)){
      ge_ranks[i, ] <- names(sort(Eg_means[i, ], decreasing = TRUE))
      dimnames(ge_ranks) <- list(rownames(Eg_means), c(1:ncol(Eg_means)))
    }

    return(
      list(
          "ge_means"  = ge_means
        , "ge_means1" = ge_means1
        , "g_means"   = g_means
        , "e_means"   = e_means
        , "ge_ranks"  = ge_ranks
      ))
  }
