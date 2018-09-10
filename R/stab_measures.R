if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
     "GE.Mean"
    , "Var"
    , "Mean"
    , "Ecov"
    , "ShuklaVar"
    ))
}


#' @name    stab_measures
#' @aliases stab_measures
#' @title  Stability Measures for Genotypes by Environment Interaction (GEI)
#' @description Stability Measures for Genotypes by Environment Interaction (GEI)
#'
#' @param .data  data.frame
#' @param .y     Response Variable
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
#' @import tibble
#' @import ggplot2
#' @import scales
#' @importFrom matrixStats rowSds rowVars
#' @importFrom reshape2 acast
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.StabMeasures <- stab_measures(
#'                 .data  = ge_data
#'                , .y    = Yield
#'                , .gen  = Gen
#'                , .env  = Env
#'                )
#' Yield.StabMeasures
#'
#'
#'
stab_measures <- function(.data, .y, .gen, .env) {
  UseMethod("stab_measures")
}

#' @export
#' @rdname stab_measures

stab_measures.default <-
  function(.data, .y, .gen, .env){

    Y   <- deparse(substitute(.y))
    G   <- deparse(substitute(.gen))
    E   <- deparse(substitute(.env))

    g <- length(levels(.data[[G]]))
    e <- length(levels(.data[[E]]))



    g_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G)) %>%
      dplyr::summarize(Mean = mean(!!rlang::sym(Y)))


    ge_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G), !!rlang::sym(E)) %>%
      dplyr::summarize(GE.Mean = mean(!!rlang::sym(Y))) %>%
      tidyr::spread(key = E, value = GE.Mean)

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

    ge_means1 <- as.matrix(ge_means[, -1])

    GE.Means    <- ge_means1
    Gen.Means   <- rowMeans(GE.Means)
    Gen.SS      <- (e-1)*matrixStats::rowVars(GE.Means)
    Gen.Var     <- matrixStats::rowVars(GE.Means)
    Gen.CV      <- matrixStats::rowSds(GE.Means)/rowMeans(GE.Means)*100
    Ecovalence  <- matrixStats::rowVars(ge_effects)*(e-1)
    GE.SS       <- sum(Ecovalence)
    GE.MSE      <- sum(Ecovalence)/((g-1)*(e-1))
    Shukla.Var  <- Ecovalence * g/((g-2)*(e-1)) - GE.MSE/(g-2)


    StabMeasures <-
      tibble::as_tibble(data.frame(
          g_means
        , GenSS     = Gen.SS
        , Var       = Gen.Var
        , CV        = Gen.CV
        , Ecov      = Ecovalence
        , ShuklaVar = Shukla.Var
        ))

    MeanVarPlot <-
      ggplot(data = StabMeasures, mapping = aes(x = Var, y = Mean)) +
      geom_point() +
      geom_text(aes(label = G), size = 2.5, vjust = 1.25, colour = "black") +
      geom_hline(yintercept = mean(StabMeasures$Mean), linetype = "dotdash") +
      labs(x = "Variance", y = "Mean") +
      scale_x_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      scale_y_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      theme_bw()

    MeanEcovPlot <-
      ggplot(data = StabMeasures, mapping = aes(x = Ecov, y = Mean)) +
      geom_point() +
      geom_text(aes(label = G), size = 2.5, vjust = 1.25, colour = "black") +
      geom_hline(yintercept = mean(StabMeasures$Mean), linetype = "dotdash") +
      labs(x = "Ecovalence", y = "Mean") +
      scale_x_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      scale_y_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      theme_bw()

    MeanShuklaVarPlot <-
      ggplot(data = StabMeasures, mapping = aes(x = ShuklaVar, y = Mean)) +
      geom_point() +
      geom_text(aes(label = G), size = 2.5, vjust = 1.25, colour = "black") +
      geom_hline(yintercept = mean(StabMeasures$Mean), linetype = "dotdash") +
      labs(x = "Shukla Variance", y = "Mean") +
      scale_x_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      scale_y_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      theme_bw()


    return(list(
            StabMeasures      = StabMeasures
          , MeanVarPlot       = MeanVarPlot
          , MeanEcovPlot      = MeanEcovPlot
          , MeanShuklaVarPlot = MeanShuklaVarPlot
          ))
  }
