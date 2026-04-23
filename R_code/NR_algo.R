
# Function init_NR for guessing initial value for NR algo based on initial parameter associated to best loglik among nrep ones 
# # INPUT: 
# nrep: number of initial parameters to test (keep it sufficiently high) (default value: 30)
# data: list of observed sequences
# # OUTPUT
# the initial parameter as (log eta, log epsilon) associated to best log-likelihood
init_NR = function(nrep = 30, data, alpha, verbose = FALSE) {
  init = matrix(NA, nrep, 2)
  loglik = rep(NA, nrep)
  for (rep in 1:nrep) {
    init[rep,] = log(runif(2, 0.01, 1.5))
    if(verbose) {cat("rep=", rep, "\n")}
    loglik[rep] = grw_der(init[rep,], data, alpha, Hess = FALSE, post.x = FALSE)$loglik
  }
  return(init[which.max(loglik),])
}

# Function optim_NR for optimisation with the Newton-Raphson algorithm 
# INPUT 
# param: inital value. Parameter in same format as for grw_der
# data: dataset  in same format as for grw_der
# niter: maximal number of iterations 
# tol: stop when max(abs(difference between new and old iteration)) drops below tol
# tol.score: tolerance on the difference between score and zero
# verbose: if TRUE, parameter estimates and log-likelihood are returned at each iteration
optim_NR = function(param, data, alpha, niter = 20, tol = 1e-6, tol.score = 1e-10, verbose = FALSE) {
  for (iter in 1:niter) {
    res = grw_der(param, data, alpha, Hess = TRUE, post.x = FALSE)
    new = param - solve(res$Hessian) %*% res$Score  
    if(verbose) {cat("param=",param,"loglik=",res$loglik, "tol=", max(abs(new-param)), "maxScore=", abs(max(res$Score)), "\n")}
    if(is.nan(res$loglik)) {stop("Convergence issue. Try more initial values for init_NR ?")}
    #check convergence
    if(max(abs(new-param)) < tol) {
      if(abs(max(res$Score)) > tol.score) {warning("Non zero score. Try more initial values for init_NR?")}
      break;
    }
    param = new
    if(iter == niter) warning("Maximal number of iterations reached. Try more initial values for init_NR?")
  }
  return(list(param=param, Score=res$Score, Hessian=res$Hessian)) 
}
