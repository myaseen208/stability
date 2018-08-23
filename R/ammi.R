#' @name    ammi
#' @aliases ammi
#' @title  Additive Main Effects and Multiplicative Interaction (AMMI)
#' @description Performs Additive Main Effects and Multiplicative Interaction (AMMI) Analysis for Genotypes by Environment Interaction (GEI)
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
#' @importFrom matrixStats rowSds rowVars
#' @importFrom stats anova as.formula ave coef confint lm pf terms aov model.tables
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.ammi <-
#'      ammi(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'       )
#' Yield.ammi
#'
#'
#'
ammi <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("ammi")
}

#' @export
#' @rdname ammi

ammi.default <-
  function(.data, .y, .rep, .gen, .env){

    Y   <- quo_name(enquo(.y))
    Rep <- quo_name(enquo(.rep))
    G   <- quo_name(enquo(.gen))
    E   <- quo_name(enquo(.env))

    g <- length(levels(.data$G))
    e <- length(levels(.data$E))
    r <- length(levels(.data$Rep))
    Min.G.E <- min(g, e)


    GE.ANOVA <-
      add_anova(
          .data = .data
        , .y    = .data$Y
        , .rep  = .data$Rep
        , .gen  = .data$G
        , .env  = .data$E
      )[[1]]


     fm2 <- aov(.data$Y ~ .data$E*.data$G + Error(.data$E/.data$Rep))
     GE.Effs <- t(model.tables(fm2, type = "effects", cterms = ".data$E:.data$G")$tables$".data$E:.data$G")

    SVD       <- svd(GE.Effs)
    D         <- diag(SVD$d[1:Min.G.E])
    GenM      <- SVD$u %*% sqrt(D)
    EnvM      <- SVD$v %*% sqrt(D)
    Ecolnumb  <- c(1:Min.G.E)
    Ecolnames <- paste0("PC", Ecolnumb)
    dimnames(GenM) <- list(levels(.data$G), Ecolnames)
    dimnames(EnvM) <- list(levels(.data$E), Ecolnames)


    SVD.Values <- SVD$d
    PC.No <- length(SVD.Values)
    PC.SS <- (SVD.Values^2)*r
    PC.Percent.SS <- PC.SS/sum(PC.SS)*100

    DF.AMMI <- rep(0, Min.G.E)
    Acum    <- DF.AMMI
    MS.AMMI <- DF.AMMI
    F.AMMI  <- DF.AMMI
    Prob.F  <- DF.AMMI
    Acumula  <- 0

    for (i in 1:Min.G.E)
    {
      DF <- (g - 1) + (e - 1) - (2 * i - 1)
      if (DF <= 0) break

      DF.AMMI[i] <- DF
      Acumula    <- Acumula + PC.Percent.SS[i]
      Acum[i]    <- Acum[i] + Acumula
      MS.AMMI[i] <- PC.SS[i]/DF.AMMI[i]
      F.AMMI[i]  <- MS.AMMI[i]/GE.ANOVA["Residual", "Mean Sq"]
      Prob.F[i]  <- 1 - pf(F.AMMI[i], DF.AMMI[i], GE.ANOVA["Residual", "Df"])
    }

    GE.Table <-
      tibble::as_tibble(data.frame(
         "PC"        = paste("PC", 1:PC.No, sep = "")
        , "SingVal"   = SVD.Values
        , "Percent"   = PC.Percent.SS
        , "accumPerc" = Acum
        , "Df"        = DF.AMMI
        , "SS"        = PC.SS
        , "Mean Sq"   = MS.AMMI
        , "F value"   = F.AMMI
        , "Pr(>F)"    = Prob.F
      ))


    GE.SS <- (t(as.vector(GE.Effs)) %*% as.vector(GE.Effs))*r

    Out <-
      list(
            anova     = GE.ANOVA
          , pc.anova  = GE.Table
      )
    return(Out)
  }
