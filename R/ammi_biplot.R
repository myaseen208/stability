if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
    "Mean"
  , "PC1"
    ))
}

#' @name    ammi_biplot
#' @aliases ammi_biplot
#' @title  Additive Main Effects and Multiplicative Interaction (AMMI) Biplot
#' @description Plots Additive Main Effects and Multiplicative Interaction (AMMI) for Genotypes by Environment Interaction (GEI)
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
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#'
#' @import rlang
#' @import ggplot2
#' @import ggfortify
#' @import scales
#' @importFrom matrixStats rowSds rowVars
#' @importFrom stats anova as.formula ave coef confint lm pf terms aov model.tables prcomp
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#'      ammi_biplot(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'       )
#'
#'
#'
ammi_biplot <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("ammi_biplot")
}
#
#' @export
#' @rdname ammi_biplot

ammi_biplot.default <-
  function(.data, .y, .rep, .gen, .env){

    Y   <- deparse(substitute(.y))
    Rep <- deparse(substitute(.rep))
    G   <- deparse(substitute(.gen))
    E   <- deparse(substitute(.env))

    g_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G)) %>%
      dplyr::summarize(Mean = mean(!!rlang::sym(Y)))

    fm2 <- aov(.data[[Y]] ~ .data[[E]]*.data[[G]] + Error(.data[[E]]/.data[[Rep]]))
    GE.Effs <- t(model.tables(fm2, type = "effects", cterms = ".data[[E]]*.data[[G]]")$tables$".data[[E]]:.data[[G]]")

      GE.AMMI <- stats::prcomp(GE.Effs, scale. = FALSE)

     aami.biplot <-
       ggplot2::autoplot(
           object         = GE.AMMI
         , label          = TRUE
        , loadings.label = TRUE
      ) +
      scale_x_continuous(sec.axis = dup_axis()) +
      scale_y_continuous(sec.axis = dup_axis()) +
      theme_bw()

    MeanPCs <-
      data.frame(
           g_means
         , GE.AMMI$x
      ) %>%
      tibble::as_tibble()


    MeanPC1Plot <-
      ggplot(data = MeanPCs, mapping = aes(x = PC1, y = Mean)) +
      geom_point() +
      geom_text(aes(label = G), size = 2.5, vjust = 1.25, colour = "black") +
      geom_vline(xintercept = 0, linetype = "dotdash") +
      geom_hline(yintercept = mean(MeanPCs$Mean), linetype = "dotdash") +
      labs(x = "PC1", y = "Mean") +
      scale_x_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      scale_y_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      theme_bw()

    return(list(

        aami.biplot = aami.biplot
        , MeanPC1Plot = MeanPC1Plot
          ))

  }
