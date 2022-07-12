#' Estimate average partial effect via augmented balancing
#' @param X the input features
#' @param Y the observed response (real valued)
#' @param W the effect variable (real valued)
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param alpha tuning parameter for glmnet
#' @param standardize whether non-binary features should be noramlized
#' @param solver convex optimizer used by CVXR for minimax weights
#' @param verbose whether the optimizer should print progress information
#' @return ATE estimate with standard error estimate. Also returns ``linear''
#'         point estimate of the form sum gamma_i Yi, as in Donoho (1994), for comparison.
#' @export average_partial_effect
average_partial_effect = function(X, Y, W,
                                  zeta = 0.5, alpha = 1, standardize = TRUE,
                                  solver = c("MOSEK", "ECOS", "SCS"),
                                  verbose = FALSE) {
    solver = match.arg(solver)
    # scale covariate matrix
    scl = apply(X, 2, sd, na.rm = TRUE)
    is.binary = apply(X, 2, \(xx) sum(xx == 0) + sum(xx == 1) == length(xx))
    scl[is.binary] = 1
    X.scl = scale(X, center = FALSE, scale = scl)

    # Compute regression adjustment
    lasso.out = rlasso2(X.scl, Y, W, alpha)
    tau.hat = lasso.out$tau_hat
    w.hat   = lasso.out$p_hat
    y.hat   = lasso.out$m_hat

    # Compute balancing weights
    gamma = balance_minimax(X.scl, W, zeta, solver = solver, verbose = verbose)

    # Compute point estimate and standard errors
    m.hat          = y.hat + (W - w.hat) * tau.hat
    point.estimate = mean(tau.hat + gamma * (Y - m.hat))
    Vhat           = mean(
                        (tau.hat - point.estimate)^2 +
                         gamma^2 * (Y - m.hat)^2
                      )
    standard.error.estimate = sqrt(Vhat / length(W))
    ret = c(
        point.estimate = point.estimate,
        standard.error.estimate = standard.error.estimate,
        linear.point.estimate = mean(gamma * Y)
    )
    return(ret)
}
