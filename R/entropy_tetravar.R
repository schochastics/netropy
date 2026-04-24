#' @title Tetravariate Entropy
#' @description Computes tetravariate entropies, expected conditional entropies,
#' and expected conditional joint entropies for all quadruples of variables in a
#' multivariate discrete data set.
#' @param dat dataframe with rows as observations and columns as variables.
#' Variables must be categorical with finite range spaces.
#' @param dec number of decimals used for rounding the entropy values.
#' Default is 2.
#' @return A dataframe with one row for each ordered decomposition of four
#' variables into predictors and conditioning variables. The columns are:
#' \item{X}{first variable in the pair of interest.}
#' \item{Y}{second variable in the pair of interest.}
#' \item{Z}{first conditioning variable.}
#' \item{U}{second conditioning variable.}
#' \item{H_XYZU}{tetravariate entropy \emph{H(X,Y,Z,U)}.}
#' \item{EH_U_XYZ}{expected conditional entropy \emph{EH(U|X,Y,Z)}.}
#' \item{EH_Z_XYU}{expected conditional entropy \emph{EH(Z|X,Y,U)}.}
#' \item{EJ_XY_ZU}{expected conditional joint entropy \emph{EJ(X,Y|Z,U)}.}
#' @details For four variables \emph{X}, \emph{Y}, \emph{Z}, and \emph{U}, the
#' tetravariate entropy is denoted \emph{H(X,Y,Z,U)}. The expected conditional
#' entropies are computed as
#' \deqn{EH(U|X,Y,Z) = H(X,Y,Z,U) - H(X,Y,Z)}
#' and
#' \deqn{EH(Z|X,Y,U) = H(X,Y,Z,U) - H(X,Y,U).}
#' The expected conditional joint entropy is computed as
#' \deqn{EJ(X,Y|Z,U) = H(X,Z,U) + H(Y,Z,U) - H(Z,U) - H(X,Y,Z,U).}
#' This quantity measures deviation from conditional independence of the form
#' \emph{X \perp Y | Z,U}. Smaller values indicate weaker conditional dependence.
#' @author Termeh Shafie
#' @seealso \code{\link{entropy_trivar}}, \code{\link{entropy_bivar}},
#' \code{\link{prediction_power}}
#' @references Frank, O., & Shafie, T. (2016). Multivariate entropy analysis of
#' network data. \emph{Bulletin of Sociological Methodology/Bulletin de
#' Méthodologie Sociologique}, 129(1), 45-63.
#' @examples
#' # use internal data set
#' data(lawdata)
#'
#' # extract node attributes
#' df_att <- lawdata[[4]]
#'
#' # data editing:
#' # 1. discretize 'years' and 'age' into three approximately balanced groups
#' # 2. recode selected variables so categories start at 0
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
#' # compute tetravariate entropy quantities for five selected variables
#' tetravariate_entropy(
#'   dat = att_var[, c("gender", "years", "age", "office", "practice")]
#' )
#'
#' @export

entropy_tetravar <- function(dat, dec = 2) {

  entropy_emp <- function(x) {
    tab <- table(x)
    p <- as.vector(tab) / sum(tab)
    p <- p[p > 0]
    sum(p * log2(1 / p))
  }

  vars <- colnames(dat)
  out <- list()
  k <- 1

  for (quad in combn(vars, 4, simplify = FALSE)) {

    xy_pairs <- combn(quad, 2, simplify = FALSE)

    for (xy in xy_pairs) {

      zu <- setdiff(quad, xy)

      x <- xy[1]
      y <- xy[2]
      z <- zu[1]
      u <- zu[2]

      h_xyzu <- entropy_emp(dat[c(x, y, z, u)])
      h_xyz  <- entropy_emp(dat[c(x, y, z)])
      h_xyu  <- entropy_emp(dat[c(x, y, u)])
      h_xzu  <- entropy_emp(dat[c(x, z, u)])
      h_yzu  <- entropy_emp(dat[c(y, z, u)])
      h_zu   <- entropy_emp(dat[c(z, u)])

      out[[k]] <- data.frame(
        X = x,
        Y = y,
        Z = z,
        U = u,
        H_XYZU = h_xyzu,
        EH_U_XYZ = h_xyzu - h_xyz,
        EH_Z_XYU = h_xyzu - h_xyu,
        EJ_XY_ZU = h_xzu + h_yzu - h_zu - h_xyzu
      )

      k <- k + 1
    }
  }

  res <- do.call(rbind, out)
  res[, 5:8] <- round(res[, 5:8], dec)
  rownames(res) <- NULL

  res
}
