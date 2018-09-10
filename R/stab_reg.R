if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(
      "Mean"
    , "GEMean"
    , "Slope"
  ))
}

#' @name    stab_reg
#' @aliases stab_reg
#' @title Individual Regression for each Genotype
#' @description Individual Regression for each Genotype in Genotypes by Environment Interaction (GEI)
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
#' @importFrom lme4 lmList
#' @importFrom magrittr %>%
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.StabReg <-
#'         stab_reg(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'           )
#'
#' Yield.StabReg
#'
stab_reg <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("stab_reg")
}

#' @export
#' @rdname stab_reg

stab_reg.default <-
  function(.data, .y, .rep, .gen, .env){

    Y   <- deparse(substitute(.y))
    Rep <- deparse(substitute(.rep))
    G   <- deparse(substitute(.gen))
    E   <- deparse(substitute(.env))

    g <- length(levels(.data[[G]]))
    e <- length(levels(.data[[E]]))
    r <- length(levels(.data[[Rep]]))


    g_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G)) %>%
      dplyr::summarize(Mean = mean(!!rlang::sym(Y)))

    DataNew  <-
      .data %>%
        dplyr::group_by(!!rlang::sym(G), !!rlang::sym(E)) %>%
        dplyr::summarize(GEMean = mean(!!rlang::sym(Y))) %>%
        dplyr::group_by(!!rlang::sym(E)) %>%
        dplyr::mutate(EnvMean = mean(GEMean))

    names(DataNew) <- c("Gen", "Env", "GEMean", "EnvMean")

    IndvReg <- lme4::lmList(GEMean ~ EnvMean|Gen, data = DataNew)
    IndvRegFit <- summary(IndvReg)

    StabIndvReg <-
      data.frame(
           g_means
        , "Slope" = coef(IndvRegFit)[ , , 2][ ,1]
        , "LCI"   = confint(IndvReg)[ , ,2][ ,1]
        , "UCI"   = confint(IndvReg)[ , ,2][ ,2]
        , "R.Sqr" = IndvRegFit$r.squared
        , "RMSE"  = IndvRegFit$sigma
        , "SSE"   = IndvRegFit$sigma^2*IndvRegFit$df[ ,2]
        , "Delta" = IndvRegFit$sigma^2*IndvRegFit$df[ ,2]/r
      ) %>%
      tibble::as_tibble()


    MeanSlopePlot <-
      ggplot(data = StabIndvReg, mapping = aes(x = Slope, y = Mean)) +
      geom_point() +
      geom_text(aes(label = G), size = 2.5, vjust = 1.25, colour = "black") +
      geom_vline(xintercept = 1, linetype = "dotdash") +
      geom_hline(yintercept = mean(StabIndvReg$Mean), linetype = "dotdash") +
      labs(x = "Slope", y = "Mean") +
      scale_x_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      scale_y_continuous(sec.axis = dup_axis(), labels = scales::comma) +
      theme_bw()


    return(
      list(
           StabIndvReg   = StabIndvReg
         , MeanSlopePlot = MeanSlopePlot
        ))
  }

