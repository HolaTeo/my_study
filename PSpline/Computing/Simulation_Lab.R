#' generate simulated response for multilevel splines
#'
#' Generates simulated response for multilevel splines
#'
#' @author YD Hwang \email{yhwang@@g.skku.edu} and ER Lee \email{erlee@@skku.edu}
#' @importFrom stats coef glm lm rbinom rnorm vcov
#' @param J  number of 'data' intervals
#' @param mod  underlying model; either `lm` or `glm`
#' @param x_sigma  design matrix sigma
#' @param e_sigma  error variance - around the mean function; data level.
#' @param z_sigma  error variance around my surface; structural level.
#' @param N_s the minimum sample size for each interval.0
#' @param N_m  the maximum sample size for each interval; default = 200.
#' @return returns a list described above.
#' @format list(x_list = x_list, y_list = y_list, e_list = e_list, true_mu = mu, z = z)
#' \describe{
#'   \item{x_list}{the length-J list of design matrices. The nrow of each element is between N_s and N_m}
#'   \item{y_list}{the length-J list of response vectors. The length of each element is between N_s and N_m.}
#'   \item{e_list}{the length-J list of error vectors. The length of each element is between N_s and N_m.}
#'   \item{true_mu}{the true mu vector of length J}
#'   \item{z}{the grid vector of length J}
#' }
#' @export

generate_response <- function(J, mod, e_sigma = 1, x_sigma = 1, z_sigma = 0.5, N_s, N_m = 200) {
  
  # currently the data interval (z interval) is set to be between -3 and 3.
  
  n <- sample(N_s:N_m, J, replace = TRUE)
  
  # smooth surface: z is the grid sequence and mu is the generated smooth function.
  z <- seq(from = -3, to = 3, length.out = J)
  mu <- z^2 - 10 * cos(2 * pi * z)  # "true" surface.
  
  beta_1 <- mu + rnorm(J, 0, z_sigma)  # slope
  beta_0 <- 0  # intercept
  
  x_list <- lapply(n, rnorm, mean = 0, sd = x_sigma)
  e_list <- lapply(n, rnorm, mean = 0, sd = e_sigma)
  
  # outcome generation function; gives 'y' list given e, beta_0, beta_1, and
  # x (design matrix)
  # for glm: logit link binary p(y = 1) = 1/(1 + exp(-beta_0 - beta_1 * x - e)
  # for lm: ordinary linear model structure y = xb + e
  if (mod == "glm") {
    y_list <- mapply(function(x, e, b, beta_0 = 0)
      rbinom(length(x), 1, 1/(1 + exp(-beta_0 - b * x - e))),
      x = x_list, e = e_list, b = beta_1)
  }
  if (mod == "lm") {
    y_list <- mapply(function(x, e, b, beta_0 = 0)
      beta_0 + b * x + e, x = x_list, e = e_list, b = beta_1)
  }
  list(x_list = x_list, y_list = y_list, e_list = e_list, true_mu = mu, z = z)
}

#' Builds ``granular'' data
#'
#' obtains the regression slope and its variance
#' certainly not optimal but this step shouldn't take long regardless
#' @author YD Hwang \email{yhwang@@g.skku.edu} and ER Lee \email{erlee@@skku.edu}
#' @param x_k design matrix
#' @param y_k response vector
#' @param mod underlying model; either `lm` or `glm`
#' @export


granular <- function(x_k, y_k, mod) {
  # summarizing the regression part
  if (mod == "glm")
    fit_lm <- glm(y_k ~ x_k, family = "binomial")
  if (mod == "lm")
    fit_lm <- lm(y_k ~ x_k)
  
  kth_beta_hat <- coef(fit_lm)[2]
  kth_var <- diag(vcov(fit_lm))[2]
  grain_out <- list(kth_beta_hat, kth_var)
  grain_out
}

#' Generates kerel matrix
#'
#' Generates kernel matrix of J by J, where J = length(z) for multilevel splines
#' certainly not optimal but this step shouldn't take long regardless.
#' Used the formulation from Reinsch (1967).
#' @author YD Hwang \email{yhwang@@g.skku.edu} and ER Lee \email{erlee@@skku.edu}
#' @param z Mid-interval value vector, it is safe to assume this to be equi-distant, but in principle it doesn't have to be. it's not tested though.
#' @export

make_K <- function(z) {
  J <- length(z)
  Del <- matrix(0, nrow = J - 2, ncol = J)
  W <- matrix(0, nrow = J - 2, ncol = J - 2)
  h <- diff(z)
  for (l in 1:(J - 2)) {
    Del[l, l] <- 1/h[l]
    Del[l, (l + 1)] <- -1/h[l] - 1/h[(l + 1)]
    Del[l, (l + 2)] <- 1/h[(l + 1)]
    W[(l - 1), l] <- W[l, (l - 1)] <- h[l]/6
    W[l, l] <- (h[l] + h[l + 1])/3
  }
  K <- t(Del) %*% solve(W) %*% Del
  K
}


#' Main EM function
#'
#' Running EM for multilevel splines
#' certainly not optimal...
#' @author YD Hwang \email{yhwang@@g.skku.edu} and ER Lee \email{erlee@@skku.edu}
#' @param beta_hat_vec data vector of length J
#' @param V covariance matrix of size J by J
#' @param K kernel matrix from `make_K`
#' @param lambda tuning parameter
#' @param maxit maximum iteration number
#' @export


main_EM <- function(beta_hat_vec, V, K, lambda, maxit = 500) {
  
  # parameter initilization
  eps <- 1000  # convergence tracker
  tol <- 1e-05  # convergence threshold
  sigma2_m <- mean(diag(V))
  J <- length(beta_hat_vec)
  mu_m <- rep(mean(beta_hat_vec), J)
  I <- diag(J)
  iter <- 1
  
  while (eps > tol & iter <= maxit) {
    # .. EM starts here
    mu_m_old <- mu_m
    sigma2_m_old <- sigma2_m  # current sigma^2
    
    Vst <- solve(solve(V) + (1/sigma2_m) * diag(J))  # Vst
    D_m <- Vst %*% solve(V)  #D_m <- part_cov %*% V
    mu_m <- solve(D_m + lambda * K) %*% D_m %*% beta_hat_vec
    
    S_lambda <- solve(I %*% D_m %*% I + lambda * K) %*% I %*% D_m
    effective_df <- sum(diag(S_lambda))
    
    sigma2_m <- mean((beta_hat_vec - mu_m)^2)
    eps <- sum(abs(mu_m - mu_m_old)) + abs(sigma2_m_old - sigma2_m)
    iter <- iter + 1
    if (iter == maxit) {
      cat("for lambda =", lambda, "max iteration reached; may need to double check \n")
    }
  }  # end of EM .. convergence reached.
  
  BIC <- sum((beta_hat_vec - mu_m)^2)/(J^(1 - effective_df/J))
  GCV <- sum((beta_hat_vec - mu_m)^2)/(J - effective_df)^2 * J
  
  EM_out <- list(mu = mu_m, S_lambda = S_lambda, sigma2 = sigma2_m, BIC = BIC, GCV = GCV)
  EM_out
}

#' Naive strawman
#'
#' Running naive splines
#' @author YD Hwang \email{yhwang@@g.skku.edu} and ER Lee \email{erlee@@skku.edu}
#' @param beta_hat_vec data vector of length J
#' @param K kernel matrix from `make_K`
#' @param lambda tuning parameter
#' @export

naive_ss <- function(beta_hat_vec, lambda, K) {
  
  J <- length(beta_hat_vec)
  I <- diag(J)
  S_lambda <- solve(I + lambda * K)
  f_hat <- S_lambda %*% beta_hat_vec
  
  eff_df <- sum(diag(S_lambda))
  
  GCV <- sum((beta_hat_vec - f_hat)^2)/(J - eff_df)^2 * J
  BIC <- log(mean((beta_hat_vec - f_hat)^2)) + eff_df * log(J)/J
  
  out <- list(mu = f_hat, S_lambda = S_lambda, BIC = BIC, GCV = GCV)
  out
}


#' Generates simulated response for multilevel splines -- test function #2
#'
#' @author YD Hwang \email{yhwang@@g.skku.edu} and ER Lee \email{erlee@@skku.edu}
#' @importFrom stats coef glm lm rbinom rnorm vcov
#' @param J  number of 'data' intervals
#' @param mod  underlying model; either `lm` or `glm`
#' @param x_sigma  design matrix sigma
#' @param e_sigma  error variance - around the mean function; data level.
#' @param z_sigma  error variance around my surface; structural level.
#' @param N_s the minimum sample size for each interval.
#' @param N_m  the maximum sample size for each interval; default = 200.
#' @return returns a list described above.
#' @format list(x_list = x_list, y_list = y_list, e_list = e_list, true_mu = mu, z = z)
#' \describe{
#' This function is supposed to be combined with the other generation function.. but later.
#'   \item{x_list}{the length-J list of design matrices. The nrow of each element is between N_s and N_m}
#'   \item{y_list}{the length-J list of response vectors. The length of each element is between N_s and N_m.}
#'   \item{e_list}{the length-J list of error vectors. The length of each element is between N_s and N_m.}
#'   \item{true_mu}{the true mu vector of length J}
#'   \item{z}{the grid vector of length J}
#' }
#' @export

generate_response_smooth <- function(J, mod, e_sigma = 1, x_sigma = 1, z_sigma = 0.5, N_s, N_m = 200) {
  
  # currently the data interval (z interval) is set to be between 0 and 1
  
  n <- sample(N_s:N_m, J, replace = TRUE)
  
  # smooth surface: z is the grid sequence and mu is the generated smooth function.
  z <- seq(from = 0, to = 1, length.out = J)
  mu <- sin(12*(z + 0.2)) / (z + 0.2)  # "true" surface.
  
  beta_1 <- mu + rnorm(J, 0, z_sigma)  # slope
  beta_0 <- 0  # intercept
  
  x_list <- lapply(n, rnorm, mean = 0, sd = x_sigma)
  e_list <- lapply(n, rnorm, mean = 0, sd = e_sigma)
  
  # outcome generation function; gives 'y' list given e, beta_0, beta_1, and
  # x (design matrix)
  # for glm: logit link binary p(y = 1) = 1/(1 + exp(-beta_0 - beta_1 * x - e)
  # for lm: ordinary linear model structure y = xb + e
  if (mod == "glm") {
    y_list <- mapply(function(x, e, b, beta_0 = 0)
      rbinom(length(x), 1, 1/(1 + exp(-beta_0 - b * x - e))),
      x = x_list, e = e_list, b = beta_1)
  }
  if (mod == "lm") {
    y_list <- mapply(function(x, e, b, beta_0 = 0)
      beta_0 + b * x + e, x = x_list, e = e_list, b = beta_1)
  }
  list(x_list = x_list, y_list = y_list, e_list = e_list, true_mu = mu, z = z)
}



library(devtools)
devtools::install_github("nclJoshCowley/PSplinesR", force = T)

#' Main EM function using pspline
#' 
#' 
main_EM_p <- function(beta_hat_vec, V, B, D, lambda, maxit = 1000) {
  
  # parameter initilization
  eps <- 1000  # convergence tracker
  tol <- 1e-05  # convergence threshold
  sigma2_m <- mean(diag(V))
  J <- length(beta_hat_vec)
  mu_m <- rep(mean(beta_hat_vec), J)
  I <- diag(J)
  iter <- 1
  
  while (eps > tol & iter <= maxit) {
    # .. EM starts here
    mu_m_old <- mu_m
    sigma2_m_old <- sigma2_m  # current sigma^2
    
    Vst <- solve(solve(V) + (1/sigma2_m) * diag(J))  # Vst
    D_m <- Vst %*% solve(V)  #D_m <- part_cov %*% V
    mu_m <- B%*%solve(t(B)%*%D_m%*%B + lambda * t(D)%*%D) %*%t(B)%*% D_m %*% beta_hat_vec
    
    S_lambda <- B%*%solve(t(B)%*%D_m%*%B + lambda * t(D)%*%D) %*%t(B)%*% D_m
    effective_df <- sum(diag(S_lambda))
    
    sigma2_m <- mean((beta_hat_vec - mu_m)^2)
    eps <- sum(abs(mu_m - mu_m_old)) + abs(sigma2_m_old - sigma2_m)
    iter <- iter + 1
    if (iter == maxit) {
      cat("for lambda =", lambda, "max iteration reached; may need to double check \n")
    }
  }  # end of EM .. convergence reached.
  
  BIC <- sum((beta_hat_vec - mu_m)^2)/(J^(1 - effective_df/J))
  GCV <- sum((beta_hat_vec - mu_m)^2)/(J - effective_df)^2 * J
  
  EM_out <- list(mu = mu_m, S_lambda = S_lambda, sigma2 = sigma2_m, BIC = BIC, GCV = GCV)
  EM_out
}


#' Naive strawman using p-spline
#' 
#' 
naive_ss_p <- function(beta_hat_vec, lambda, B, D) {
  
  J <- length(beta_hat_vec)
  I <- diag(J)
  S_lambda <- B%*%solve(t(B)%*%B + lambda * t(D)%*%D) %*%t(B)
  f_hat <- S_lambda %*% beta_hat_vec
  
  eff_df <- sum(diag(S_lambda))
  
  GCV <- sum((beta_hat_vec - f_hat)^2)/(J - eff_df)^2 * J
  BIC <- log(mean((beta_hat_vec - f_hat)^2)) + eff_df * log(J)/J
  
  out <- list(mu = f_hat, S_lambda = S_lambda, BIC = BIC, GCV = GCV)
  out
}


# Simulation results



## Ex1_lm_single

#### Example 1

### LM

## Single

RMSE_SINGLE_LM_1 <- NULL
SD_SINGLE_LM_1 <- NULL

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_single_lm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      RMSE_single_lm_1[j] <- sqrt((1/50)*t(y$true_mu-unlist(betahat[,1]))%*%(y$true_mu-unlist(betahat[,1])))
      
    }
    
    RMSE_SINGLE_lm_1 <- mean(RMSE_single_lm_1)
    SD_SINGLE_lm_1 <- sd(RMSE_single_lm_1)
    RMSE_SINGLE_LM_1 <- cbind(RMSE_SINGLE_LM_1,RMSE_SINGLE_lm_1)
    SD_SINGLE_LM_1 <- cbind(SD_SINGLE_LM_1,SD_SINGLE_lm_1)
    
  }
  
}

RMSE_SINGLE_LM_1 ; SD_SINGLE_LM_1

## Ex1_lm_naive


#### Example 1

### LM

## Naive

RMSE_NAIVE_LM_1 <- NULL
SD_NAIVE_LM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_naive_lm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_lm_1[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVE_lm_1 <- mean(RMSE_naive_lm_1)
    SD_NAIVE_lm_1 <- sd(RMSE_naive_lm_1)
    RMSE_NAIVE_LM_1 <- cbind(RMSE_NAIVE_LM_1,RMSE_NAIVE_lm_1)
    SD_NAIVE_LM_1 <- cbind(SD_NAIVE_LM_1,SD_NAIVE_lm_1)
    
  }
}

RMSE_NAIVE_LM_1 ; SD_NAIVE_LM_1



## Ex1_lm_naive(pspline)



#### Example 1

### LM

## Naive

RMSE_NAIVEP_LM_1 <- NULL
SD_NAIVEP_LM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_naive_lm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_lm_1[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVEP_lm_1 <- mean(RMSE_naive_lm_1)
    SD_NAIVEP_lm_1 <- sd(RMSE_naive_lm_1)
    RMSE_NAIVEP_LM_1 <- cbind(RMSE_NAIVEP_LM_1,RMSE_NAIVEP_lm_1)
    SD_NAIVEP_LM_1 <- cbind(SD_NAIVE_LM_1,SD_NAIVE_lm_1)
    
  }
}
RMSE_NAIVEP_LM_1 ; SD_NAIVEP_LM_1



## Ex1_lm_multilevel



#### Example 1

### LM

## Multilevel

RMSE_MULTI_LM_1 <- NULL
SD_MULTI_LM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_multi_lm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_lm_1[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTI_lm_1 <- mean(RMSE_multi_lm_1)
    SD_MULTI_lm_1 <- sd(RMSE_multi_lm_1)
    RMSE_MULTI_LM_1 <- cbind(RMSE_MULTI_LM_1,RMSE_MULTI_lm_1)
    SD_MULTI_LM_1 <- cbind(SD_MULTI_LM_1,SD_MULTI_lm_1)
    
  }
}

RMSE_MULTI_LM_1 ; SD_MULTI_LM_1



## Ex1_lm_multilevel(pspline)



#### E#### Example 1

### LM

## Multilevel

RMSE_MULTIP_LM_1 <- NULL
SD_MULTIP_LM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01)


for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    
    betahat <- NULL
    RMSE_multi_lm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_lm_1[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTIP_lm_1 <- mean(RMSE_multi_lm_1)
    SD_MULTIP_lm_1 <- sd(RMSE_multi_lm_1)
    RMSE_MULTIP_LM_1 <- cbind(RMSE_MULTIP_LM_1,RMSE_MULTIP_lm_1)
    SD_MULTIP_LM_1 <- cbind(SD_MULTIP_LM_1,SD_MULTIP_lm_1)
    
    
  }
}

RMSE_MULTIP_LM_1 ; SD_MULTIP_LM_1


## Ex1_glm_single



#### Example 1

### GLM

## Single

RMSE_SINGLE_GLM_1 <- NULL
SD_SINGLE_GLM_1 <- NULL

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    
    betahat <- NULL
    RMSE_single_glm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      RMSE_single_glm_1[j] <- sqrt((1/50)*t(y$true_mu-unlist(betahat[,1]))%*%(y$true_mu-unlist(betahat[,1])))
      
    }
    
    RMSE_SINGLE_glm_1 <- mean(RMSE_single_glm_1)
    SD_SINGLE_glm_1 <- sd(RMSE_single_glm_1)
    RMSE_SINGLE_GLM_1 <- cbind(RMSE_SINGLE_GLM_1,RMSE_SINGLE_glm_1)
    SD_SINGLE_GLM_1 <- cbind(SD_SINGLE_GLM_1,SD_SINGLE_glm_1)
    
    
  }
}

RMSE_SINGLE_GLM_1 ; SD_SINGLE_GLM_1




## Ex1_glm_naive



#### Example 1

### GLM

## Naive

RMSE_NAIVE_GLM_1 <- NULL
SD_NAIVE_GLM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)


for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    
    betahat <- NULL
    RMSE_naive_glm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_glm_1[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVE_glm_1 <- mean(RMSE_naive_glm_1)
    SD_NAIVE_glm_1 <- sd(RMSE_naive_glm_1)
    RMSE_NAIVE_GLM_1 <- cbind(RMSE_NAIVE_GLM_1,RMSE_NAIVE_glm_1)
    SD_NAIVE_GLM_1 <- cbind(SD_NAIVE_GLM_1,SD_NAIVE_glm_1)
    
  }
}

RMSE_NAIVE_GLM_1 ; SD_NAIVE_GLM_1 




## Ex1_glm_naive(pspline)



#### Example 1

### GLM

## Naive

RMSE_NAIVEP_GLM_1 <- NULL
SD_NAIVEP_GLM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_naive_glm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_glm_1[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVEP_glm_1 <- mean(RMSE_naive_glm_1)
    SD_NAIVEP_glm_1 <- sd(RMSE_naive_glm_1)
    RMSE_NAIVEP_GLM_1 <- cbind(RMSE_NAIVEP_GLM_1,RMSE_NAIVEP_glm_1)
    SD_NAIVEP_GLM_1 <- cbind(SD_NAIVEP_GLM_1,SD_NAIVEP_glm_1)
    
  }
}

RMSE_NAIVEP_GLM_1 ; SD_NAIVEP_GLM_1 




## Ex1_glm_multilevel



#### Example 1

### GLM

## Multilevel

RMSE_MULTI_GLM_1 <- NULL
SD_MULTI_GLM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_multi_glm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_glm_1[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTI_glm_1 <- mean(RMSE_multi_glm_1)
    SD_MULTI_glm_1 <- sd(RMSE_multi_glm_1)
    RMSE_MULTI_GLM_1 <- cbind(RMSE_MULTI_GLM_1,RMSE_MULTI_glm_1)
    SD_MULTI_GLM_1 <- cbind(SD_MULTI_GLM_1,SD_MULTI_glm_1)
    
    
  }
}

RMSE_MULTI_GLM_1 ; SD_MULTI_GLM_1 




## Ex1_glm_multilevel(pspline)



#### Example 1

### GLM

## Multilevel

RMSE_MULTIP_GLM_1 <- NULL
SD_MULTIP_GLM_1 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_multi_glm_1 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response_smooth(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_glm_1[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTIP_glm_1 <- mean(RMSE_multi_glm_1)
    SD_MULTIP_glm_1 <- sd(RMSE_multi_glm_1)
    RMSE_MULTIP_GLM_1 <- cbind(RMSE_MULTIP_GLM_1,RMSE_MULTIP_glm_1)
    SD_MULTIP_GLM_1 <- cbind(SD_MULTIP_GLM_1,SD_MULTIP_glm_1)
    
  }
}
RMSE_MULTIP_GLM_1 ; SD_MULTIP_GLM_1




## Ex2_lm_single


#### Example 2

### LM

## Single

RMSE_SINGLE_LM_2 <- NULL
SD_SINGLE_LM_2 <- NULL

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_single_lm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      RMSE_single_lm_2[j] <- sqrt((1/50)*t(y$true_mu-unlist(betahat[,1]))%*%(y$true_mu-unlist(betahat[,1])))
      
    }
    
    RMSE_SINGLE_lm_2 <- mean(RMSE_single_lm_2)
    SD_SINGLE_lm_2 <- sd(RMSE_single_lm_2)
    RMSE_SINGLE_LM_2 <- cbind(RMSE_SINGLE_LM_2,RMSE_SINGLE_lm_2)
    SD_SINGLE_LM_2 <- cbind(SD_SINGLE_LM_2,SD_SINGLE_lm_2)
    
  }
}

RMSE_SINGLE_LM_2 ; SD_SINGLE_LM_2 



## Ex2_lm_naive



#### Example 2

### LM

## Naive

RMSE_NAIVE_LM_2 <- NULL
SD_NAIVE_LM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_naive_lm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_lm_2[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVE_lm_2 <- mean(RMSE_naive_lm_2)
    SD_NAIVE_lm_2 <- sd(RMSE_naive_lm_2)
    RMSE_NAIVE_LM_2 <- cbind(RMSE_NAIVE_LM_2,RMSE_NAIVE_lm_2)
    SD_NAIVE_LM_2 <- cbind(SD_NAIVE_LM_2,SD_NAIVE_lm_2) 
    
  }
}

RMSE_NAIVE_LM_2 ; SD_NAIVE_LM_2




## Ex2_lm_naive(pspline)



#### Example 2

### LM

## Naive

RMSE_NAIVEP_LM_2 <- NULL
SD_NAIVEP_LM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)


for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_naive_lm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_lm_2[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVEP_lm_2 <- mean(RMSE_naive_lm_2)
    SD_NAIVEP_lm_2 <- sd(RMSE_naive_lm_2)
    RMSE_NAIVEP_LM_2 <- cbind(RMSE_NAIVEP_LM_2,RMSE_NAIVEP_lm_2)
    SD_NAIVEP_LM_2 <- cbind(SD_NAIVEP_LM_2,SD_NAIVEP_lm_2)
    
  }
}
RMSE_NAIVEP_LM_2 ; SD_NAIVEP_LM_2



## Ex2_lm_multilevel



#### Example 2

### LM

## Multilevel

RMSE_MULTI_LM_2 <- NULL
SD_MULTI_LM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.0005)


for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_multi_lm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_lm_2[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTI_lm_2 <- mean(RMSE_multi_lm_2)
    SD_MULTI_lm_2 <- sd(RMSE_multi_lm_2)
    RMSE_MULTI_LM_2 <- cbind(RMSE_MULTI_LM_2,RMSE_MULTI_lm_2)
    SD_MULTI_LM_2 <- cbind(SD_MULTI_LM_2,SD_MULTI_lm_2)
    
  }
}

RMSE_MULTI_LM_2 ; SD_MULTI_LM_2 



## Ex2_lm_multilevel(pspline)



#### Example 2

### LM

## Multilevel

RMSE_MULTIP_LM_2 <- NULL
SD_MULTIP_LM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01)


for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_multi_lm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "lm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "lm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_lm_2[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTIP_lm_2 <- mean(RMSE_multi_lm_2)
    SD_MULTIP_lm_2 <- sd(RMSE_multi_lm_2)
    RMSE_MULTIP_LM_2 <- cbind(RMSE_MULTIP_LM_2,RMSE_MULTIP_lm_2)
    SD_MULTIP_LM_2 <- cbind(SD_MULTIP_LM_2,SD_MULTIP_lm_2)
    
    
  }
}

RMSE_MULTIP_LM_2 ; SD_MULTIP_LM_2 





## Ex2_glm_single



#### Example 2

### GLM

## Single

RMSE_SINGLE_GLM_2 <- NULL
SD_SINGLE_GLM_2 <- NULL

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_single_glm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      RMSE_single_glm_2[j] <- sqrt((1/50)*t(y$true_mu-unlist(betahat[,1]))%*%(y$true_mu-unlist(betahat[,1])))
      
    }
    
    RMSE_SINGLE_glm_2 <- mean(RMSE_single_glm_2)
    SD_SINGLE_glm_2 <- sd(RMSE_single_glm_2)
    RMSE_SINGLE_GLM_2 <- cbind(RMSE_SINGLE_GLM_2,RMSE_SINGLE_glm_2)
    SD_SINGLE_GLM_2 <- cbind(SD_SINGLE_GLM_2,SD_SINGLE_glm_2)
    
  }
}

RMSE_SINGLE_GLM_2 ; SD_SINGLE_GLM_2



## Ex2_glm_naive



#### Example 2

### GLM

## Naive

RMSE_NAIVE_GLM_2 <- NULL
SD_NAIVE_GLM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)


for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_naive_glm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss(beta_hat_vec = unlist(betahat[,1]), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_glm_2[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVE_glm_2 <- mean(RMSE_naive_glm_2)
    SD_NAIVE_glm_2 <- sd(RMSE_naive_glm_2)
    RMSE_NAIVE_GLM_2 <- cbind(RMSE_NAIVE_GLM_2,RMSE_NAIVE_glm_2)
    SD_NAIVE_GLM_2 <- cbind(SD_NAIVE_GLM_2,SD_NAIVE_glm_2)
    
  }
}

RMSE_NAIVE_GLM_2 ; SD_NAIVE_GLM_2



## Ex2_glm_naive(pspline)



#### Example 2

### GLM

## Naive

RMSE_NAIVEP_GLM_2 <- NULL
SD_NAIVEP_GLM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001,  0.000001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_naive_glm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,naive_out$GCV)
      }
      
      naive_out <- naive_ss_p(beta_hat_vec = unlist(betahat[,1]), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_naive_glm_2[j] <- sqrt((1/50)*t(y$true_mu-naive_out$mu)%*%(y$true_mu-naive_out$mu))
      
    }
    
    RMSE_NAIVEP_glm_2 <- mean(RMSE_naive_glm_2)
    SD_NAIVEP_glm_2 <- sd(RMSE_naive_glm_2)
    RMSE_NAIVEP_GLM_2 <- cbind(RMSE_NAIVEP_GLM_2,RMSE_NAIVEP_glm_2)
    SD_NAIVEP_GLM_2 <- cbind(SD_NAIVEP_GLM_2,SD_NAIVEP_glm_2)
    
    
  }
}

RMSE_NAIVEP_GLM_2 ; SD_NAIVEP_GLM_2



## Ex2_glm_multilevel



#### Example 2

### GLM

## Multilevel

RMSE_MULTI_GLM_2 <- NULL
SD_MULTI_GLM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001)

for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_multi_glm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      K <- make_K(y$z)
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), K = K, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_glm_2[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTI_glm_2 <- mean(RMSE_multi_glm_2)
    SD_MULTI_glm_2 <- sd(RMSE_multi_glm_2)
    RMSE_MULTI_GLM_2 <- cbind(RMSE_MULTI_GLM_2,RMSE_MULTI_glm_2)
    SD_MULTI_GLM_2 <- cbind(SD_MULTI_GLM_2,SD_MULTI_glm_2)
    
  }
}


RMSE_MULTI_GLM_2 ; SD_MULTI_GLM_2 



## Ex2_glm_multilevel(pspline)



#### Example 2

### GLM

## Multilevel

RMSE_MULTIP_GLM_2 <- NULL
SD_MULTIP_GLM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01)


for(t in c(2,4,8)){
  
  for(n in c(50,100)){
    
    betahat <- NULL
    RMSE_multi_glm_2 <- NULL
    
    for(j in 1:200){
      
      betahat <- NULL
      y <-  generate_response(J = 50, mod = "glm", e_sigma = t, N_s = n)
      
      for(i in 1:50){
        
        results <- granular(unlist(y$x_list[[i]]), unlist(y$y_list[[i]]), mod = "glm")
        betahat <- rbind(betahat,results)
        
      }
      
      B <- PSplinesR::GetBSpline(x = y$z, deg = 3, IntKnots = y$z[-c(1, length(y$z))], ExtKnots = y$z[c(1, length(y$z))] )
      D <- PSplinesR::GetDiffMatrix(dim(B)[1], 2) 
      
      GCV_vec <- NULL
      for(k in 1:length(lambda)){
        EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[k])
        GCV_vec <- rbind(GCV_vec,EM_out$GCV)
      }
      
      EM_out <- main_EM_p(beta_hat_vec = unlist(betahat[,1]), V = diag(unlist(betahat[,2])), B = B, D = D, lambda = lambda[which.min(GCV_vec)])
      
      RMSE_multi_glm_2[j] <- sqrt((1/50)*t(y$true_mu-EM_out$mu)%*%(y$true_mu-EM_out$mu))
      
    }
    
    RMSE_MULTIP_glm_2 <- mean(RMSE_multi_glm_2)
    SD_MULTIP_glm_2 <- sd(RMSE_multi_glm_2)
    RMSE_MULTIP_GLM_2 <- cbind(RMSE_MULTIP_GLM_2,RMSE_MULTIP_glm_2)
    SD_MULTIP_GLM_2 <- cbind(SD_MULTIP_GLM_2,SD_MULTIP_glm_2)
    
  }
}

RMSE_MULTIP_GLM_2 ; SD_MULTIP_GLM_2



## Result

RMSE_SINGLE_LM_1 ; SD_SINGLE_LM_1
RMSE_NAIVE_LM_1 ; SD_NAIVE_LM_1
RMSE_NAIVEP_LM_1 ; SD_NAIVEP_LM_1
RMSE_MULTI_LM_1 ; SD_MULTI_LM_1
RMSE_MULTIP_LM_1 ; SD_MULTIP_LM_1

RMSE_SINGLE_GLM_1 ; SD_SINGLE_GLM_1
RMSE_NAIVE_GLM_1 ; SD_NAIVE_GLM_1 
RMSE_NAIVEP_GLM_1 ; SD_NAIVEP_GLM_1
RMSE_MULTI_GLM_1 ; SD_MULTI_GLM_1
RMSE_MULTIP_GLM_1 ; SD_MULTIP_GLM_1

RMSE_SINGLE_LM_2 ; SD_SINGLE_LM_2 
RMSE_NAIVE_LM_2 ; SD_NAIVE_LM_2
RMSE_NAIVEP_LM_2 ; SD_NAIVEP_LM_2
RMSE_MULTI_LM_2 ; SD_MULTI_LM_2
RMSE_MULTIP_LM_2 ; SD_MULTIP_LM_2

RMSE_SINGLE_GLM_2 ; SD_SINGLE_GLM_2
RMSE_NAIVE_GLM_2 ; SD_NAIVE_GLM_2
RMSE_NAIVEP_GLM_2 ; SD_NAIVEP_GLM_2
RMSE_MULTI_GLM_2 ; SD_MULTI_GLM_2
RMSE_MULTIP_GLM_2 ; SD_MULTIP_GLM_2


