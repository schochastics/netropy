#' @title Divergence Tests of Goodness of Fit
#' @description Performs divergence-based goodness-of-fit tests for discrete data,
#' including tests of uniformity, pairwise independence, conditional independence,
#' and nested model comparisons.
#' @param dat dataframe with rows as observations and columns as variables.
#' Variables must be categorical with finite range spaces.
#' @param var character name of a variable in \code{dat} to test for uniformity.
#' @param var1 character name of first variable.
#' @param var2 character name of second variable.
#' @param var_cond optional character vector of conditioning variables.
#' @param model_full list containing \code{D} and \code{df} for the full model.
#' @param model_reduced list containing \code{D} and \code{df} for the reduced model.
#' @param alpha significance level (default 0.05).
#' @param dec number of decimals for rounding (default 3).
#' @param use_approx_cv logical; if TRUE uses approximation for critical value,
#' otherwise uses exact chi-square quantile.
#' @return Dataframe with divergence \emph{D}, chi-square statistic, degrees of freedom,
#' critical value, and decision.
#' @details
#' The function implements four types of tests:
#'
#' \strong{1. Uniformity}
#' \deqn{D = \log r_X - H(X)}
#'
#' \strong{2. Pairwise Independence}
#' \deqn{D = H(X) + H(Y) - H(X,Y)}
#'
#' \strong{3. Conditional Independence}
#' \deqn{D = H(X,Z) + H(Y,Z) - H(Z) - H(X,Y,Z)}
#'
#' \strong{4. Nested Model Comparison}
#' \deqn{2n[D(p,q_2) - D(p,q_1)]}
#'
#' The test statistic is given by
#' \deqn{2nD}
#' and is approximately chi-square distributed.
#'
#' Smaller divergence values indicate better model fit.
#' @author Termeh Shafie
#' @seealso \code{\link{joint_entropy}}, \code{\link{entropy_trivar}}
#' @references Frank, O., & Shafie, T. (2016). Multivariate entropy analysis of network data.
#' @examples
#' # use internal dataset
#' data(lawdata)
#'
#' # extract attributes
#' df_att <- lawdata[[4]]
#'
#' # discretize variables
#' att_var <- data.frame(
#'   status    = df_att$status - 1,
#'   gender    = df_att$gender,
#'   office    = df_att$office - 1,
#'   years     = ifelse(df_att$years <= 3, 0,
#'                 ifelse(df_att$years <= 13, 1, 2)),
#'   age       = ifelse(df_att$age <= 35, 0,
#'                 ifelse(df_att$age <= 45, 1, 2)),
#'   practice  = df_att$practice,
#'   lawschool = df_att$lawschool - 1
#' )
#'
#' ## 1. Test uniformity
#' div_gof(att_var, var = "gender")
#'
#' ## 2. Test independence
#' div_gof(att_var, var1 = "status", var2 = "gender")
#'
#' ## 3. Test conditional independence
#' div_gof(att_var,
#'         var1 = "status",
#'         var2 = "gender",
#'         var_cond = "years")
#'
#' ## 4. Nested model comparison
#'  ## Compare a reduced independence model against the saturated empirical model
#' m_full <- list(D = 0, df = 0)
#'
#' m_reduced <- div_gof(att_var,
#'                     var1 = "status",
#'                     var2 = "gender")
#'
#' div_gof(att_var,
#'         model_full = m_full,
#'         model_reduced = list(D = m_reduced$D, df = m_reduced$df))
#'
#'  ## Nested comparison for conditional independence
#' m_full <- list(D = 0, df = 0)
#'
#' m_reduced <- div_gof(att_var,
#'                     var1 = "status",
#'                     var2 = "gender",
#'                     var_cond = "years")
#'
#' div_gof(att_var,
#'         model_full = m_full,
#'         model_reduced = list(D = m_reduced$D, df = m_reduced$df))
#' @export
#'

div_gof <- function(dat,
                    var = NULL,
                    var1 = NULL,
                    var2 = NULL,
                    var_cond = NULL,
                    model_full = NULL,
                    model_reduced = NULL,
                    alpha = 0.05,
                    dec = 3,
                    use_approx_cv = TRUE) {

  entropy_set <- function(dat, vars) {
    tab <- table(dat[, vars, drop = FALSE])
    p <- as.vector(tab) / sum(tab)
    p <- p[p > 0]
    sum(p * log2(1 / p))
  }

  n_levels <- function(x) length(unique(x))

  make_output <- function(test, D, df_chi2) {
    chi2_stat <- 2 * nrow(dat) * D * log(2)

    if (use_approx_cv) {
      crit_val <- df_chi2 + sqrt(8 * df_chi2)
    } else {
      crit_val <- qchisq(1 - alpha, df = df_chi2)
    }

    decision <- if (chi2_stat > crit_val) {
      "reject"
    } else {
      "cannot reject"
    }

    data.frame(
      test = test,
      D = round(D, dec),
      chi2 = round(chi2_stat, dec),
      df = df_chi2,
      critical_value = round(crit_val, dec),
      decision = decision
    )
  }

  # 1. uniformity: X ~ uniform
  if (!is.null(var) && is.null(var1) && is.null(var2)) {
    r_x <- n_levels(dat[[var]])

    D <- log2(r_x) - entropy_set(dat, var)
    df_chi2 <- r_x - 1

    return(make_output(
      test = paste0("uniformity: ", var),
      D = D,
      df_chi2 = df_chi2
    ))
  }

  # 2. pairwise independence: X independent of Y
  if (!is.null(var1) && !is.null(var2) && is.null(var_cond)) {
    r_x <- n_levels(dat[[var1]])
    r_y <- n_levels(dat[[var2]])

    D <- entropy_set(dat, var1) +
      entropy_set(dat, var2) -
      entropy_set(dat, c(var1, var2))

    df_chi2 <- (r_x - 1) * (r_y - 1)

    return(make_output(
      test = paste0(var1, " independent of ", var2),
      D = D,
      df_chi2 = df_chi2
    ))
  }

  # 3. conditional independence: X independent of Y given Z
  if (!is.null(var1) && !is.null(var2) && !is.null(var_cond)) {
    r_x <- n_levels(dat[[var1]])
    r_y <- n_levels(dat[[var2]])
    r_cond <- prod(sapply(dat[var_cond], n_levels))

    D <- entropy_set(dat, c(var1, var_cond)) +
      entropy_set(dat, c(var2, var_cond)) -
      entropy_set(dat, var_cond) -
      entropy_set(dat, c(var1, var2, var_cond))

    df_chi2 <- (r_x - 1) * (r_y - 1) * r_cond

    return(make_output(
      test = paste0(var1, " independent of ", var2,
                    " given ", paste(var_cond, collapse = ", ")),
      D = D,
      df_chi2 = df_chi2
    ))
  }

  # 4. nested model comparison
  if (!is.null(model_full) && !is.null(model_reduced)) {
    D <- model_reduced$D - model_full$D
    df_chi2 <- model_reduced$df - model_full$df

    if (df_chi2 <= 0) {
      stop(
        "Invalid nested model comparison: `model_reduced` must impose more constraints ",
        "than `model_full`, so its degrees of freedom must be larger."
      )
    }

    if (D < 0) {
      stop(
        "Invalid nested model comparison: divergence for the reduced model ",
        "should be at least as large as for the full model."
      )
    }

    return(make_output(
      test = "nested model comparison",
      D = D,
      df_chi2 = df_chi2
    ))
  }
}
