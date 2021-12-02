

## X is array, Z is matrix, theta is matrix

## Note that X's row is parameter dimension, column is the space dimension




############################################################# step 5  update v,sig2,rho

update_v = function(Z, v, Tstar, av, bv){
  T = nrow(Z)
  
  eta = rbeta(1, v+1, T);
  c = av+Tstar
  d = bv-log(eta)
  
  
  ind = sample(x = c(0, 1), size = 1, prob = c(c-1, T*d))
  
  if (ind == 0){
    v_new = rgamma(n = 1, shape = c, rate = d)
  } else {
    v_new = rgamma(n = 1, shape = c-1, rate = d)
  }
  
  v_new
}



update_sig2 = function(asig, bsig, D, W, rho, Theta){
  
  Theta_star = unique.matrix(Theta)
  
  Tstar = nrow(Theta_star)
  n = ncol(Theta)
  
  H_inv = D-rho*W
  H_inv = as.matrix(H_inv)
  
  sum_thetaH = 0
  for (j in 1:Tstar)
  {
    sum_thetaH = sum_thetaH + t(Theta_star[j,])%*%H_inv%*%as.matrix(Theta_star[j,])
  }
  
  c = asig + Tstar*n/2
  d = bsig + 0.5*sum_thetaH
  
  sig2 = rinvgamma(1,c,d)
  
  sig2
}



update_rho = function(D, W, Theta, sig2){
  
  Theta_star = unique.matrix(Theta)
  
  Tstar = nrow(Theta_star)
  n = ncol(Theta)
  
  
  myden = function(rho)
  {
    H = solve(D - rho*W)
    H_inv = D - rho*W
    H_inv = as.matrix(H_inv)
    
    sum_thetaH = 0
    for (j in 1:Tstar)
    {
      sum_thetaH = sum_thetaH + t(Theta_star[j,])%*%H_inv%*%as.matrix(Theta_star[j,])
    }
    
    return(-Tstar/2*as.numeric(determinant(H)$modulus) - 1/(2*sig2)*sum_thetaH)
  }
  
  indmeden = function(rho) (rho>0)*(rho<1)
  
  
  rho_new = arms(0.5, myden, indmeden, n.sample = 1)
  
  rho_new
  
}







############################################################# step 4  update Beta, tau2

update_Beta = function(Betam, Sigma_m, tau2, X, Y, Theta, R){
  
  
  T = nrow(Y)
  sum_xx = sum_xytheta = 0
  for (t in 1:T)
  {
    sum_xx = sum_xx + R[,,t]%*%t(R[,,t])
    sum_xytheta = sum_xytheta + R[,,t]%*%as.matrix(Y[t,] - X[t]*Theta[t,])
  }
  
  Sigma_Beta_tilde = solve(solve(Sigma_m) + (1/tau2)*sum_xx)
  Beta_tilde = Sigma_Beta_tilde%*%(solve(Sigma_m)%*%Betam + (1/tau2)*sum_xytheta)
  
  
  c(rmvnorm(1, Beta_tilde, Sigma_Beta_tilde))
  
  
  
}




update_tau2 = function(X, Y, Theta, atau, btau, Beta, R){
  
  T = nrow(Y)
  n = ncol(Y)
  sum_yxbetatheta = 0
  for (t in 1:T)
  {
    sum_yxbetatheta = sum_yxbetatheta + t(Y[t,] - X[t]*Theta[t,] - t(R[,,t])%*%Beta)%*%(Y[t,] - X[t]*Theta[t,] - t(R[,,t])%*%Beta)
  }
  
  c = atau + n*T/2
  d = btau + 0.5*sum_yxbetatheta
  
  tau2 = rinvgamma(1, c, d)
  
  
  tau2
}







############################################################# step 2 & 3  update  Theta
############################################################# used "rgsdp.cpp" to speed up 




############################################################# step 1 generate Y from Z

update_Y = function(Z, X, tau2, Theta, Beta, R){
  
  ###### our model indicating the ith and jth element of mu are independent 
  T = nrow(Z)
  n = ncol(Z)
  MU = matrix(0, nrow = T, ncol = n)
  for (t in 1:T)  MU[t,] = X[t]*Theta[t,] + t(R[,,t])%*%Beta
  
  sd = sqrt(tau2)
  
  Y = matrix(0, nrow = T, ncol = n)
  for (t in 1:T){
    
    for (j in 1:n){
      
      if (Z[t,j] == 0) {
        
        phi = pnorm(0, mean = MU[t,j], sd = sd)
        
        u = runif(1, 0, phi)
        
        Y[t,j] = qnorm(u, mean = MU[t,j], sd = sd)
        
      } else {
        
        phi_lower = pnorm(Z[t,j] - 1, mean = MU[t,j], sd = sd, log.p = T, lower.tail = TRUE)
        phi_upper = pnorm(Z[t,j], mean = MU[t,j], sd = sd, log.p = T, lower.tail = TRUE)
        
        u = runif(1, phi_lower, phi_upper)
        
        Y[t,j] = qnorm(u, mean = MU[t,j], sd = sd, log.p = T, lower.tail = TRUE)
        
        if (is.infinite(Y[t,j])){
          phi_lower = pnorm(Z[t,j] - 1, mean = MU[t,j], sd = sd, log.p = T, lower.tail = FALSE)
          phi_upper = pnorm(Z[t,j], mean = MU[t,j], sd = sd, log.p = T, lower.tail = FALSE)
          
          u = runif(1, phi_upper, phi_lower)
          
          Y[t,j] = qnorm(u, mean = MU[t,j], sd = sd, log.p = T, lower.tail = FALSE)
        }
        
        
      }
      
    }
    
  }
  
  Y
  
}



