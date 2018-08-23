#' @name    add_anova
#' @aliases add_anova
#' @title Additive ANOVA for Genotypes by Environment Interaction (GEI) model
#' @description Additive ANOVA for Genotypes by Environment Interaction (GEI) model
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
#' @importFrom stats anova as.formula ave coef confint formula lm pf terms
#'
#' @export
#'
#' @examples
#' data(ge_data)
#' YieldANOVA <-
#'      add_anova(
#'             .data = ge_data
#'           , .y    = Yield
#'           , .rep  = Rep
#'           , .gen  = Gen
#'           , .env  = Env
#'       )
#' YieldANOVA
#'
add_anova <- function(.data, .y, .rep, .gen, .env) {
  UseMethod("add_anova")
}

#' @export
#' @rdname add_anova

add_anova.default <-
  function(.data, .y, .rep, .gen, .env){
    Y   <- quo_name(enquo(.y))
    Rep <- quo_name(enquo(.rep))
    G   <- quo_name(enquo(.gen))
    E   <- quo_name(enquo(.env))

    fm1 <- formula(paste0(Y, "~", paste(E, paste(Rep, E, sep = ":"),
                                        G, paste(G, E, sep = ":"), sep = "+")))

    fm2 <- lm(terms(fm1, keep.order = TRUE), data = .data)
    fm2ANOVA <- anova(fm2)
    rownames(fm2ANOVA) <- c("Env", "Rep(Env)", "Gen", "Gen:Env",
                            "Residuals")
    fm2ANOVA[1, 4] <- fm2ANOVA[1, 3]/fm2ANOVA[2, 3]
    fm2ANOVA[2, 4] <- NA
    fm2ANOVA[1, 5] <- 1 - pf(as.numeric(fm2ANOVA[1, 4]), fm2ANOVA[1,
                                                                  1], fm2ANOVA[2, 1])
    fm2ANOVA[2, 5] <- 1 - pf(as.numeric(fm2ANOVA[2, 4]), fm2ANOVA[2,
                                                                  1], fm2ANOVA[5, 1])
    class(fm2ANOVA) <- c("anova", "data.frame")
    return(list(anova = fm2ANOVA))

  }

