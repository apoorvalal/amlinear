#' Estimate minimax balancing weights
#'
#' @param X the input features
#' @param W the effect variable (real valued)
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param solver convex optimzer used by CVXR
#' @param verbose whether convex optimizer should print verbose output
#'
#' @return balancing weights
#'
#' @export balance_minimax
balance_minimax = function(X, W, zeta, solver = c("MOSEK", "ECOS", "SCS"), verbose = TRUE) {
  solver = match.arg(solver)
  n = nrow(X); p = ncol(X)
  γ = CVXR::Variable(n + 2)
  objective = ((1 - zeta) * sum(γ[1:n]^2) +         # l2 norm of weights
                     zeta * sum(γ[n + 1:2]^2)       # agg imbalance terms
              )
  constraints = list(
      sum(γ[1:n]) == 0,
       t(X) %*% γ[1:n] <= γ[n + 1],   # control mean
      -t(X) %*% γ[1:n] <= γ[n + 1],   # treatment mean
       sum(W * γ[1:n]) == 1,
       t(X) %*% (W * γ[1:n]) <=   colMeans(X) + γ[n + 2],
      -t(X) %*% (W * γ[1:n]) <= - colMeans(X) + γ[n + 2]
  )
  cvx.problem = CVXR::Problem(CVXR::Minimize(objective), constraints)
  cvx.output = solve(cvx.problem, solver = solver, verbose = verbose)
  result = cvx.output$getValue(γ)
  gamma = n * result[1:n]
  gamma
}
