#' @name    gge_biplot
#' @aliases gge_biplot
#' @title  Genotype plus Genotypes by Environment (GGE) Interaction Biplot
#' @description Plots Genotype plus Genotypes by Environment (GGE) Interaction Biplot for Genotypes by Environment Interaction (GEI)
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .rep   Replication Factor
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#'
#' @return Stability Measures
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
#'
#' @import rlang
#' @import ggplot2
#' @import ggfortify
#' @importFrom matrixStats rowSds rowVars
#' @importFrom stats anova as.formula ave coef confint lm pf terms aov model.tables prcomp
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#'      gge_biplot(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'       )
#'
#'
#'
gge_biplot <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("gge_biplot")
}
#
#' @export
#' @rdname gge_biplot

gge_biplot.default <-
  function(.data, .y, .rep, .gen, .env){

    Y   <- deparse(substitute(.y))
    Rep <- deparse(substitute(.rep))
    G   <- deparse(substitute(.gen))
    E   <- deparse(substitute(.env))

    ge_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G), !!rlang::sym(E)) %>%
      dplyr::summarize(GE.Mean = mean(!!rlang::sym(Y))) %>%
      tidyr::spread(key = !!rlang::sym(E), value = GE.Mean)

    ge_means1 <- as.matrix(ge_means[, -1])
    rownames(ge_means1) <- c(ge_means[, 1])[[1]]

    gge_effects <-
      sweep(
        x      = ge_means1
        , MARGIN = 2
        , STATS  = colMeans(ge_means1)
      )

    GGE.AMMI <- stats::prcomp(gge_effects, scale. = FALSE)

    gge.biplot <-
      ggplot2::autoplot(
          object         = GGE.AMMI
        , label          = TRUE
        , loadings.label = TRUE
      ) +
      scale_x_continuous(sec.axis = dup_axis()) +
      scale_y_continuous(sec.axis = dup_axis()) +
      theme_bw()

    return(list(
          gge.biplot = gge.biplot
          ))

  }
