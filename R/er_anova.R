#' @name    er_anova
#' @aliases er_anova
#' @title  Eberhart & Russel’s Model ANOVA
#' @description ANOVA of Eberhart & Russel’s Model
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
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#' @import rlang
#' @import tibble
#' @importFrom lme4 lmList
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.er_anova <-
#'          er_anova(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'           )
#' Yield.er_anova
#'
#'
er_anova <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("er_anova")
}

#' @export
#' @rdname er_anova

er_anova.default <-
  function(.data, .y, .rep, .gen, .env){

    Y   <- enquo(.y)
    Rep <- enquo(.rep)
    G   <- enquo(.gen)
    E   <- enquo(.env)

    g <- length(levels(.data$G))
    e <- length(levels(.data$E))
    r <- length(levels(.data$Rep))

    # Individual Regression
    StabIndvReg <-
            stab_reg(
                .data = .data
              , .y    = !! Y
              , .rep  = !! Rep
              , .gen  = !! G
              , .env  = !! E
              )$StabIndvReg


fm1 <- lm(
            formula = terms(.data$Y ~ .data$E + .data$Rep:.data$E + .data$G + .data$G:.data$E, keep.order = TRUE)
          , data = .data
          )
fm1ANOVA <- anova(fm1)
rownames(fm1ANOVA) <- c("Env", "Rep(Env)", "Gen", "Gen:Env", "Residuals")
fm1ANOVA[1, 4] <- fm1ANOVA[1, 3]/fm1ANOVA[2, 3]
fm1ANOVA[2, 4] <- fm1ANOVA[2, 3]/fm1ANOVA[5, 3]
fm1ANOVA[1, 5] <- 1-pf(as.numeric(fm1ANOVA[1, 4]), fm1ANOVA[1, 1], fm1ANOVA[2, 1])
fm1ANOVA[2, 5] <- 1-pf(as.numeric(fm1ANOVA[2, 4]), fm1ANOVA[2, 1], fm1ANOVA[5, 1])
class(fm1ANOVA) <- c("anova", "data.frame")

    Df <-
      c(
          g*e-1
        , g-1
        , g*(e-1)
        , 1
        , g-1
        , g*(e-2)
        , rep(x = e-2, times = g)
        , e*g*(r- 1)
      )

    PooledError <- fm1ANOVA["Residuals", "Sum Sq"]/r
    TotalSS <- (fm1ANOVA["Env", "Sum Sq"] + fm1ANOVA["Rep(Env)", "Sum Sq"] + fm1ANOVA["Gen:Env", "Sum Sq"])/r
    GenSS   <- fm1ANOVA["Rep(Env)", "Sum Sq"]/r
    EnvGESS <- (fm1ANOVA["Env", "Sum Sq"] + fm1ANOVA["Gen:Env", "Sum Sq"])/r
    EnvL    <- fm1ANOVA[1, "Sum Sq"]/r
    GELSS   <- EnvGESS - EnvL - sum(StabIndvReg$Delta)
    PooledDev  <- sum(StabIndvReg$Delta)
    GensSS  <- StabIndvReg$Delta


    SS <-
      c(
          TotalSS
        , GenSS
        , EnvGESS
        , EnvL
        , GELSS
        , PooledDev
        , GensSS
        , PooledError
      )

    MS <- SS/Df

    F <-
      c(
        NA
        , MS[2]/MS[6]
        , NA
        , NA
        , MS[5]/MS[6]
        , NA
        , MS[7:(length(MS) - 1)]/MS[length(MS)]
        , NA
      )

    PLines <- 1 - pf(F[7:(length(MS) - 1)], Df[7], Df[length(Df)])
    pval <- c(NA, 1 - pf(F[2], Df[2], Df[6]), NA, NA, 1 - pf(F[5], Df[5], Df[6]), NA, PLines, NA)

    ANOVA <- data.frame(Df, `Sum Sq` = SS, `Mean Sq` = MS,
                        `F value` = F, `Pr(>F)` = pval, check.names = FALSE)

    rownames(ANOVA) <-
      c(
        "Total"
        , "Gen"
        , "Env + (Gen x Env)"
        , "  Env (linear)"
        , "  Gen x Env(linear)"
        , "  Pooled deviation"
        , paste0("    ", levels(.data$G))
        , "Pooled error"
      )

    class(ANOVA) <- c("anova", "data.frame")

    return(list(
      fm1ANOVA
      , er_anova    = ANOVA
    )
    )
  }

