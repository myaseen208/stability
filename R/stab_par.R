#' @name    stab_par
#' @aliases stab_par
#' @title  Stability Parameters for Genotypes by Environment Interaction (GEI)
#' @description Stability Parameters for Genotypes by Environment Interaction (GEI)
#'
#' @param .data   data.frame
#' @param .y      Response Variable
#' @param .rep    Replication Factor
#' @param .gen    Genotypes Factor
#' @param .env    Environment Factor
#' @param alpha   Level of Significance, default is 0.1
#' @param .envCov Environmental Covariate, default is NULL
#'
#' @return Stability Parameters
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
#' @importFrom stats anova as.formula ave coef confint lm pf terms qf qt
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.StabPar <-
#'    stab_par(
#'             .data   = ge_data
#'           , .y      = Yield
#'           , .rep    = Rep
#'           , .gen    = Gen
#'           , .env    = Env
#'           , alpha   = 0.1
#'           , .envCov = NULL
#' )
#'
#' Yield.StabPar
#'
stab_par <- function(.data, .y, .rep, .gen, .env, alpha = 0.1, .envCov = NULL) {
  UseMethod("stab_par")
}

#' @export
#' @rdname stab_par


stab_par.default <-
  function(.data, .y, .rep, .gen, .env, alpha = 0.1, .envCov = NULL) {

    Y   <- deparse(substitute(.y))
    Rep <- deparse(substitute(.rep))
    G   <- deparse(substitute(.gen))
    E   <- deparse(substitute(.env))
    EnvCov   <- deparse(substitute(.envCov))

    g <- length(levels(.data[[G]]))
    e <- length(levels(.data[[E]]))
    r <- length(levels(.data[[Rep]]))

    fm1 <- lm(formula = terms(.data[[Y]] ~ .data[[E]] + .data[[Rep]]:.data[[E]] +
                                .data[[G]] + .data[[G]]:.data[[E]], keep.order = TRUE))
    fm1ANOVA <- anova(fm1)
    rownames(fm1ANOVA) <- c("Env", "Rep(Env)", "Gen", "Gen:Env",
                            "Residuals")
    fm1ANOVA[1, 4] <- fm1ANOVA[1, 3]/fm1ANOVA[2, 3]
    fm1ANOVA[2, 4] <- NA
    fm1ANOVA[1, 5] <- 1 - pf(as.numeric(fm1ANOVA[1, 4]), fm1ANOVA[1,
                                                                  1], fm1ANOVA[2, 1])
    fm1ANOVA[2, 5] <- 1 - pf(as.numeric(fm1ANOVA[2, 4]), fm1ANOVA[2,
                                                                  1], fm1ANOVA[5, 1])
    class(fm1ANOVA) <- c("anova", "data.frame")

    MSE <- fm1ANOVA["Residuals", "Mean Sq"]

    g_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G)) %>%
      dplyr::summarize(Mean = mean(!!rlang::sym(Y)))

    ge_means0 <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G), !!rlang::sym(E)) %>%
      dplyr::summarize(GE.Mean = mean(!!rlang::sym(Y))) %>%
      tidyr::spread(key = E, value = GE.Mean)

    ge_means1 <- as.matrix(ge_means0[, -1])
    rownames(ge_means1) <- c(ge_means0[, 1])[[1]]


    x <- ge_means1

     if(is.null(.data[[EnvCov]]) == TRUE) {
       EnvCov <- colMeans(x) - mean(x)
       }
     else
       {
         EnvCov <- EnvCov
       }


    ge_means <-
      .data %>%
      dplyr::group_by(!!rlang::sym(G), !!rlang::sym(E)) %>%
      dplyr::summarize(GE.Mean = mean(!!rlang::sym(Y)))
    names(ge_means) <- c("G", "E", "GE.Mean")

    fm2ANOVA <- anova(lm(formula = GE.Mean ~ G*E, data = ge_means))

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


    Gen.Mean <- rowMeans(x)
    GenSS    <- rowVars(x)*(e - 1)
    Gen.Var  <- rowVars(x)
    Gen.CV   <- rowSds(x)/rowMeans(x)*100
    Ecov     <- rowVars(ge_effects)*(e - 1)
    GE.SS    <- sum(Ecov)
    GE.MSE   <- sum(Ecov)/((g-1)*(e-1))

    Sigma   <- g/((g-2)*(e-1)) * Ecov - GE.MSE/(g-2)
    Beta    <- ge_effects %*% EnvCov/sum(EnvCov^2)
    New.ge_effects <-
          ge_effects -
         matrix(data = EnvCov, nrow=g, ncol=e, byrow=TRUE) *
         matrix(data = Beta, nrow = g, ncol = e, byrow = FALSE)

    SP      <- (g/((g-2)*(e-2)))*(rowSums(New.ge_effects^2)- (sum(rowSums(New.ge_effects^2))/(g*(g-1))))
    SSRes  <- ((g - 1)*(e - 2)*sum(SP))/g
    SSHetro   <- fm2ANOVA["G:E", "Sum Sq"] - SSRes


     Df <- c(sum(fm2ANOVA[-4, "Df"]), fm2ANOVA[-4, "Df"], fm2ANOVA["G", "Df"], fm2ANOVA["G:E", "Df"] - fm2ANOVA["G", "Df"], fm2ANOVA["G", "Df"]^2*(fm2ANOVA["E", "Df"] + 1))
     SumS  <- c(sum(fm2ANOVA[-4, "Sum Sq"]), fm2ANOVA[-4, "Sum Sq"], SSHetro, SSRes, NA)*r
     MeanSS <- c(NA, SumS[2:6]/Df[2:6], MSE)

    F <-
      c(
        NA
        , MeanSS[2]/MeanSS[6]
        , MeanSS[3]/MeanSS[7]
        , MeanSS[4]/MeanSS[7]
        , MeanSS[5]/MeanSS[6]
        , MeanSS[6]/MeanSS[7]
        , NA
      )

    pval <-
      c(
        NA
        , 1 - pf(F[2], Df[2], Df[6])
        , 1 - pf(F[3], Df[3], Df[7])
        , 1 - pf(F[4], Df[4], Df[7])
        , 1 - pf(F[5], Df[5], Df[6])
        , 1 - pf(F[6], Df[6], Df[7])
        , NA
      )

    ANOVA <- data.frame(Df, `Sum Sq` = SumS, `Mean Sq` = MeanSS,
                        `F value` = F, `Pr(>F)` = pval, check.names = FALSE)

    rownames(ANOVA) <-
      c(
        "Total"
        , "Gen"
        , "Env"
        , "Gen x Env"
        , "  Heterogeneity"
        , "  Residual"
        , "Pooled error"
      )

    class(ANOVA) <- c("anova", "data.frame")



    dfShukla <- e*(g - 1)*(r - 1) # df for Shukla Test
    FM0 <- qf(1 - alpha, e-1, dfShukla)
    F05 <- qf(0.95, e-1, dfShukla)
    F01 <- qf(0.99, e-1, dfShukla)
    MS    <- SSRes/(fm2ANOVA["Gen", "Df"]^2*(fm2ANOVA["Env", "Df"]+1))
    DMV05 <- qt(0.95, dfShukla) * sqrt(2 * MSE/(r * e))
    MES   <-  mean(x)
    FS    <- Sigma*r/MSE
    FSS   <- SP*r/MSE

    NN <-
      ifelse(
        test = FS < F05
        , yes  = "ns"
        , no   =  ifelse(
          test = FS < F01 & FS >= F05
          , yes  = "*"
          , no = "**"
        )
      )


    MMM <-
      ifelse(
        test = FSS < F05
        , yes  = "ns"
        , no   =  ifelse(
          test = FSS < F01 & FSS >= F05
          , yes  = "*"
          , no = "**"
        )
      )

    Stb.Results <-
      tibble::as_tibble(data.frame(
         g_means
        , GenSS = GenSS/r
        , Var   = Gen.Var/r
        , CV    = Gen.CV/r
        , Ecov  = Ecov/r
        , GE.SS = GE.SS/r
        , GE.MSE = GE.MSE/r
        , "Sigma" = Sigma/r
        , "."     = NN
        , "SP"    = SP/r
        , ".."     = MMM
      ))

    # Simultaneous Selection
    F1 <- as.numeric(as.character(
      cut(
        x = FS
        , breaks = c(-Inf, FM0, F05, F01,  Inf)
        , labels = c(0, 2, -4, -8)
        , right=TRUE
      )
    ))


    MV1 <- as.numeric(as.character(
      cut(
        x = Gen.Mean
        , breaks = c(-Inf, MES-2*DMV05, MES-DMV05, MES, MES+DMV05, MES+2*DMV05, Inf)
        , labels = c(-3, -2, -1, 1, 2, 3)
        , right=FALSE
      )
    ))

    R    <- rank(x = Gen.Mean, na.last = TRUE, ties.method = "min")
    GY   <- R + MV1
    GYS  <- GY + F1
    SGYS <- sum(GYS)
    MGYS <- SGYS/g
    GYY  <- ifelse(test = GYS > MGYS, yes = "+", no =  "-")

    SimulSel <-
      tibble::as_tibble(data.frame(
        g_means, R, MV1, GY, Sigma/r, F1, GYS, GYY
        ))

    names(SimulSel) <- c("Gen", "Mean", "Rank", "Adjustment", "Adj.Rank", "Sigma", "Stab.Rating", "YSi", "Select")

    return(list(
          ANOVA     = ANOVA
        , StabPar   = Stb.Results
        , SimultSel = SimulSel
    ))
  }
