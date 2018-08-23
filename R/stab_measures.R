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

    Y   <- enquo(.y)
    G   <- enquo(.gen)
    E   <- enquo(.env)

    g <- length(levels(.data$G))
    e <- length(levels(.data$E))
    r <- length(levels(.data$Rep))

    g_means <-
          ge_means(
                .data  = .data
               , .y    = !! Y
               , .gen  = !! G
               , .env  = !! E
               )$g_means

    names(g_means) <- c("G", "Mean")

     ge_means <-
       ge_means(
         .data  = .data
         , .y    = !! Y
         , .gen  = !! G
         , .env  = !! E
       )$ge_means


    ge_means1 <-
      ge_means(
        .data  = .data
        , .y    = !! Y
        , .gen  = !! G
        , .env  = !! E
      )$ge_means1

    gge_effects <-
      ge_effects(
        .data  = .data
        , .y    = !! Y
        , .gen  = !! G
        , .env  = !! E
      )$gge_effects


    ge_effects <-
      ge_effects(
        .data  = .data
        , .y    = !! Y
        , .gen  = !! G
        , .env  = !! E
      )$ge_effects

    GE.Means    <- ge_means1
    Gen.Means   <- rowMeans(GE.Means)
    Gen.SS      <- (e-1)*rowVars(GE.Means)
    Gen.Var     <- rowVars(GE.Means)
    Gen.CV      <- rowSds(GE.Means)/rowMeans(GE.Means)*100
    Ecovalence  <- rowVars(ge_effects)*(e-1)
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
