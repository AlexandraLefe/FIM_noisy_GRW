# This file contains the main function, grw_der, which computed the log-likelihood of a noisy GRW, its gradient and its Hessian matrix

# INPUT
#' @param param: vector of size 2. Parameter as (log eta, log epsilon)
#' @param data: list of observed sequences in the same format as output z of function "simu" 
#' @param Hess: boolean. If TRUE, the Hessian of the log-likelihood is returned (default value: TRUE)
#' @param post.x: boolean. If TRUE, the list of estimated mean and variance of latent states are returned (default value = FALSE)

# OUTPUT: list of 
#' @return loglik: the log-likelihood of the model parameter
#' @return grad: the score, ie. the gradient of the log-likelihood
#' @return Hess (if HESS=TRUE as input): the Hessian matrix of the log-likelihood
#' @return post.x (if post.x=TRUE as input): the list of 
# M: a list of same length as data. Each slot contains the vector of mean of the latent states
# S2: a list of same length as data. Each slot contains the vector of variance of the latent states
grw_der = function(param, data, alpha, Hess = TRUE, post.x = FALSE) {
  
  # initialization of the results
  res.loglik = 0
  res.grad = rep(0,length(param))
  if(Hess) H11 = H21 = H22 = 0.0  
  if(post.x) {x = list(M = vector("list", length(data)), S2 = vector("list", length(data)))}
  # if(post.x) {hat.x = list(M = matrix(NA,m,n), S2 = matrix(NA,m,n))}
  
  eta = exp(param[1])
  epsilon = exp(param[2])
  
  aux.e = exp(-2*param) # stored as it will be often used
  
  # forward backward and their gradient
  for (j in 1:length(data)) { # sequence number j 
    z = data[[j]] 
    n = length(z)
    
    K = L = mu = sigma = nu = tau = rep(NA,n)
    A = B = array(NA, dim = c(n,2,3))
    
    # # forward 
    K[1] = dnorm(z[1], alpha, sqrt(epsilon^2+eta^2),log=TRUE)
    mu[1] = (z[1]*epsilon^2 + alpha*eta^2) / (epsilon^2 + eta^2)
    sigma[1] = sqrt((epsilon^2*eta^2) / (epsilon^2+eta^2))
    
    # careful product '*' differently handled by arrays and matrices. Safer to distribute by hand
    tmp = rbind(c(z[1]^2, -2*z[1], 1) * aux.e[1], 
                c(alpha^2, -2*alpha, 1) * aux.e[2])  
    tmp[,1] = tmp[,1] - 1
    A[1,,] = tmp
    
    for (i in 2:n) {
      K[i] = K[i-1] + dnorm(z[i], mu[i-1], sqrt(epsilon^2 + eta^2 + sigma[i-1]^2), log=TRUE)
      mu[i] = ((mu[i-1]) * eta^2 + z[i] * (sigma[i-1]^2 + epsilon^2)) / (eta^2 + sigma[i-1]^2 + epsilon^2)
      sigma[i] = sqrt((eta^2 * (sigma[i-1]^2 + epsilon^2)) / (eta^2 + sigma[i-1]^2 + epsilon^2))
    }
    # P, Q, R2 stored as they will be later used 
    # Lower complexity with vectorial computation and new loop on i versus all within same loop i 
    P = sigma[-n]^2 / (epsilon^2 + sigma[-n]^2)
    Q = mu[-n] * epsilon^2 / (epsilon^2 + sigma[-n]^2)
    R2 = epsilon^2 * sigma[-n]^2 / (epsilon^2 + sigma[-n]^2)
    # New recursion to compute A
    for (i in 2:n) {
      tmp = rbind(c(z[i]^2, -2*z[i], 1) * aux.e[1], 
                  c(Q[i-1]^2+R2[i-1], 2*(P[i-1]-1)*Q[i-1], (P[i-1]-1)^2) * aux.e[2])
      tmp[,1] = tmp[,1] - 1
      # TODO OPTIMZE product and sum ! Sparse matrices (package if it exists or by hand)
      A[i,,] = tmp + A[i-1,,] %*% rbind(c(1,0,0), 
                                        c(Q[i-1],P[i-1],0), 
                                        c(Q[i-1]^2+R2[i-1], 2*P[i-1]*Q[i-1], P[i-1]^2)) 
    }
    
    # log-likelihood and its gradient 
    res.loglik = res.loglik + K[n] 
    grad = A[n,,] %*% c(1, mu[n], mu[n]^2+sigma[n]^2)
    res.grad = res.grad + grad
    
    if(!Hess & !post.x) next # next sequence
    
    # backward 
    L[n] = -dnorm(0, 0, 1e5, log = TRUE)
    nu[n] = 0
    tau[n] = 1e5
    B[n,,] = matrix(0,2,3)
    
    for (i in n:2) {
      L[i-1] = L[i] + dnorm(z[i],nu[i],sqrt(eta^2+tau[i]^2),log=TRUE)
      nu[i-1] = (nu[i]*eta^2 + z[i]*tau[i]^2) / (eta^2+tau[i]^2)
      tau[i-1] = sqrt(epsilon^2+eta^2*tau[i]^2/(eta^2+tau[i]^2))
      omega2 = eta^2*tau[i]^2 / (eta^2+tau[i]^2)
      P.tilde = omega2 / (epsilon^2 + omega2)
      Q.tilde = nu[i-1] * epsilon^2 / (epsilon^2 + omega2)
      R2.tilde = P.tilde * epsilon^2
      
      tmp = rbind(c((Q.tilde-z[i])^2 + R2.tilde, 2*P.tilde*(Q.tilde-z[i]), P.tilde^2) * aux.e[1],
                  c(Q.tilde^2+R2.tilde,  2*(P.tilde-1)*Q.tilde, (P.tilde-1)^2) * aux.e[2])
      tmp[,1] = tmp[,1] - 1
      # TODO OPTIMZE: sparse matrices
      B[i-1,,] = tmp + B[i,,] %*% rbind(c(1,0,0), 
                                        c(Q.tilde,P.tilde,0), 
                                        c(Q.tilde^2+R2.tilde, 2*P.tilde*Q.tilde, P.tilde^2))
    }
    # # verification loglik using each separator and clique
    # K[n] # last FW
    # head(K + L + dnorm(mu, nu, sqrt(sigma^2 + tau^2), log=TRUE)) # separators
    # head( # cliques
    #   K[-n] + L[-1] +
    #     dnorm(z[-1], nu[-1], sqrt(eta^2 + tau[-1]^2), log=TRUE) +
    #     dnorm((z[-1] * tau[-1]^2 + nu[-1] * eta^2) / (tau[-1]^2 + eta^2), mu[-n],
    #           sqrt(eta^2 * tau[-1]^2 / (eta^2 + tau[-1]^2) + sigma[-n]^2 + epsilon^2), log=TRUE)
    # )
    
    # Vector of Mean and Variance of posterior state densities
    M = (nu*sigma^2+mu*tau^2)/(sigma^2+tau^2) 
    S2 = sigma^2*tau^2/(sigma^2+tau^2)
    if(post.x) {x$M[[j]] = M; x$S2[[j]] = S2}
    
    if(!Hess) next # next sequence or complete computation of the Hessian
    
    # Hessian
    for (i in 1:n) {
      # Matrix A+B minus the score in the first column (related to constant). Stored for later use
      aux.abg = A[i,,] + B[i,,] 
      aux.abg[,1] = aux.abg[,1] - grad
      
      # Coef H.11 (first one of the diagonal)
      tmp = c(z[i]^2, -2*z[i], 1) * aux.e[1]
      tmp[1] = tmp[1]-1
      # TODO OPTIMZE: sparse matrices
      # NB: If C created with Matrix package, total complexity is worse (empirically tested)  
      C = cbind(c(aux.abg[1,],0,0),
                c(0, aux.abg[1,],0), 
                c(0,0,aux.abg[1,])) %*% tmp - 2*aux.e[1] * c(z[i]^2, -2*z[i], 1, 0, 0)
      # library("Matrix") # for sparse matrices
      # C = Matrix(c(aux.abg[1,],0,0, 0, aux.abg[1,],0, 0,0,aux.abg[1,]), nrow=5, sparse = TRUE) %*% tmp - 2*exp(-2*param[1]) * c(z[i]^2, -2*z[i], 1, 0, 0)
      H11 = H11 + sum(t(C) %*% rbind(1, 
                                     M[i], 
                                     M[i]^2+S2[i], 
                                     M[i]^3 + 3*M[i]*S2[i], M[i]^4 + 6*M[i]^2*S2[i] + 3*S2[i]^2))
      
      # Coef H.21 (off diagonal)
      # TODO OPTIMZE: sparse matrices
      D = cbind(c(aux.abg[2,],0,0), 
                c(0, aux.abg[2,],0), 
                c(0,0,aux.abg[2,])) %*% tmp
      H21 = H21 + sum(t(D) %*% rbind(1, 
                                     M[i], 
                                     M[i]^2+S2[i], 
                                     M[i]^3 + 3*M[i]*S2[i], 
                                     M[i]^4 + 6*M[i]^2*S2[i] + 3*S2[i]^2))
    } # close loop on i for H11 and H21
    
    # Coef H.22 (second one of the diagonal)
    # # start with H22.1 (related to i=1)
    aux.abg = (A[1,,] + B[1,,])[2,]
    aux.abg[1] = aux.abg[1] - grad[2]
    tmp = c(alpha^2, -2*alpha, 1) * aux.e[2]
    tmp[1] = tmp[1]-1 
    # TODO OPTIMZE : sparse matrix
    E1 = cbind(c(aux.abg,0,0), 
               c(0, aux.abg,0), 
               c(0,0,aux.abg)) %*% tmp - 2*aux.e[2] * c(alpha^2, -2*alpha, 1, 0, 0)
    H22 = H22 + sum(t(E1) %*% rbind(1, 
                                    M[1], 
                                    M[1]^2+S2[1], 
                                    M[1]^3 + 3*M[1]*S2[1], 
                                    M[1]^4 + 6*M[1]^2*S2[1] + 3*S2[1]^2))
    
    # # pursue with H22.2 (related to i=2,..n)
    omega2 = eta^2 * tau[-1]^2 / (eta^2 + tau[-1]^2) # re-compute omega2
    N = (nu[-n] * (epsilon^2 + sigma[-n]^2) + mu[-n] * omega2) / (omega2 + epsilon^2 + sigma[-n]^2)
    T2 = omega2 * (epsilon^2 + sigma[-n]^2) / (omega2 + epsilon^2 + sigma[-n]^2)
    
    # moments 
    Moments = rbind(c(N),
                    c(N^2 + T2),
                    c(N^3 + 3*N*T2),
                    c(N^4 + 6*N^2*T2 + 3*T2^2))
    
    # # NB: second dimension of matrices A and B are indices from 0 in the paper -> add 1 
    # # NB: elements are flushed out the sum whenever possible to avoid multiple single operations
    # # store elements that will be used several times
    e00 = - A[-n,2,1] - B[-1,2,1]  + 1 + grad[2]
    e00p3e = (e00 + 3) * aux.e[2] 
    
    # part related to constant x0y0
    H22 = H22 + sum(e00) 
    # part related to x
    H22 = H22 + sum(-A[-n,2,2] * (P*N + Q))
    # part related to y
    H22 = H22 + sum(-B[-1,2,2] * N)
    # part related to x2
    H22 = H22 + sum(- (A[-n,2,3] + e00p3e) * (P^2*Moments[2,] + 2*P*Q*N + Q^2+R2))
    # part related to y2
    H22 = H22 + sum(- (B[-1,2,3] + e00p3e) * Moments[2,])
    # part related to xy
    H22 = H22 + sum(2*e00p3e * (P*Moments[2,] + Q*N))
    # part related to x3
    H22 = H22 + sum(A[-n,2,2] *
                      (P^3*Moments[3,] + 3*P^2*Q*Moments[2,] + 
                         3*P*(Q^2 + R2)*N + Q^3 + 3*Q*R2)) * aux.e[2]
    # part related to y3
    H22 = H22 + sum(B[-1,2,2] * Moments[3,]) * aux.e[2]
    # part related to x2y
    H22 = H22 + sum((-2*A[-n,2,2] + B[-1,2,2]) * 
                      (P^2*Moments[3,] + 2*P*Q*Moments[2,] + (Q^2 + R2)*N)) * aux.e[2]
    # part related to xy2
    H22 = H22 + sum((-2*B[-1,2,2] + A[-n,2,2]) *
                      (P*Moments[3,] + Q*Moments[2,])) * aux.e[2]
    # part related to x4
    H22 = H22 + sum((A[-n,2,3] + aux.e[2]) * 
                      (P^4*Moments[4,] + 4*P^3*Q*Moments[3,] + 6*P^2*(Q^2 + R2)*Moments[2,] + 
                         (4*P*Q^3 + 12*P*Q*R2)*N + Q^4 + 6*R2*Q^2 + 3*R2^2)) * aux.e[2]
    # part related to y4
    H22 = H22 + sum((B[-1,2,3] + aux.e[2]) * Moments[4,]) * aux.e[2]
    # part related to x3y
    H22 = H22 - sum((A[-n,2,3] + 2*aux.e[2]) * 
                      (P^3*Moments[4,] + 3*P^2*Q*Moments[3,] + 3*P*(Q^2 + R2)*Moments[2,] + 
                         (Q^3 + 3*Q*R2)*N)) * (2*aux.e[2])
    # part related to xy3
    H22 = H22 - sum((B[-1,2,3] + 2*aux.e[2]) * (P*Moments[4,] + Q*Moments[3,])) * (2*aux.e[2])
    # part related to x2y2
    H22 = H22 + sum((A[-n,2,3] + B[-1,2,3] + 6*aux.e[2]) * 
                      (P^2*Moments[4,] + 2*P*Q*Moments[3,] + (Q^2 + R2)*Moments[2,])) * aux.e[2]
  }
  
  res = list(loglik = res.loglik, Score = c(res.grad)) # grad = c() for converting array to vector
  if(Hess) res$Hessian = matrix(c(H11,H21,H21,H22), nrow=2)
  if(post.x) res$post.x = x
  return(res)
}
