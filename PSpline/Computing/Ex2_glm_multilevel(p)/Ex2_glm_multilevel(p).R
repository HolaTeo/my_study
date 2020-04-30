#### Example 2

### GLM

## Multilevel

RMSE_MULTI_GLM_2 <- NULL
SD_MULTI_GLM_2 <- NULL
lambda <- c(1000, 100, 10, 1, 0.1, 0.01)

# (2,50)

betahat <- NULL
RMSE_multi_glm_2 <- NULL

for(j in 1:200){
  
  betahat <- NULL
  y <-  generate_response(J = 50, mod = "glm", e_sigma = 2, N_s = 50)
  
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

RMSE_MULTI_glm_2 <- mean(RMSE_multi_glm_2)
SD_MULTI_glm_2 <- sd(RMSE_multi_glm_2)
RMSE_MULTI_GLM_2 <- cbind(RMSE_MULTI_GLM_2,RMSE_MULTI_glm_2)
SD_MULTI_GLM_2 <- cbind(SD_MULTI_GLM_2,SD_MULTI_glm_2)

# (2,100)

betahat <- NULL
RMSE_multi_glm_2 <- NULL

for(j in 1:200){
  
  betahat <- NULL
  y <-  generate_response(J = 50, mod = "glm", e_sigma = 2, N_s = 100)
  
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

RMSE_MULTI_glm_2 <- mean(RMSE_multi_glm_2)
SD_MULTI_glm_2 <- sd(RMSE_multi_glm_2)
RMSE_MULTI_GLM_2 <- cbind(RMSE_MULTI_GLM_2,RMSE_MULTI_glm_2)
SD_MULTI_GLM_2 <- cbind(SD_MULTI_GLM_2,SD_MULTI_glm_2)

# (4,50)

betahat <- NULL
RMSE_multi_glm_2 <- NULL

for(j in 1:200){
  
  betahat <- NULL
  y <-  generate_response(J = 50, mod = "glm", e_sigma = 4, N_s = 50)
  
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

RMSE_MULTI_glm_2 <- mean(RMSE_multi_glm_2)
SD_MULTI_glm_2 <- sd(RMSE_multi_glm_2)
RMSE_MULTI_GLM_2 <- cbind(RMSE_MULTI_GLM_2,RMSE_MULTI_glm_2)
SD_MULTI_GLM_2 <- cbind(SD_MULTI_GLM_2,SD_MULTI_glm_2)

# (4,100)

betahat <- NULL
RMSE_multi_glm_2 <- NULL

for(j in 1:200){
  
  betahat <- NULL
  y <-  generate_response(J = 50, mod = "glm", e_sigma = 4, N_s = 100)
  
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

RMSE_MULTI_glm_2 <- mean(RMSE_multi_glm_2)
SD_MULTI_glm_2 <- sd(RMSE_multi_glm_2)
RMSE_MULTI_GLM_2 <- cbind(RMSE_MULTI_GLM_2,RMSE_MULTI_glm_2)
SD_MULTI_GLM_2 <- cbind(SD_MULTI_GLM_2,SD_MULTI_glm_2)

# (8,50)

betahat <- NULL
RMSE_multi_glm_2 <- NULL

for(j in 1:200){
  
  betahat <- NULL
  y <-  generate_response(J = 50, mod = "glm", e_sigma = 8, N_s = 50)
  
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

RMSE_MULTI_glm_2 <- mean(RMSE_multi_glm_2)
SD_MULTI_glm_2 <- sd(RMSE_multi_glm_2)
RMSE_MULTI_GLM_2 <- cbind(RMSE_MULTI_GLM_2,RMSE_MULTI_glm_2)
SD_MULTI_GLM_2 <- cbind(SD_MULTI_GLM_2,SD_MULTI_glm_2)

# (8,100)

betahat <- NULL
RMSE_multi_glm_2 <- NULL

for(j in 1:200){
  
  betahat <- NULL
  y <-  generate_response(J = 50, mod = "glm", e_sigma = 8, N_s = 100)
  
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

RMSE_MULTI_glm_2 <- mean(RMSE_multi_glm_2)
SD_MULTI_glm_2 <- sd(RMSE_multi_glm_2)
RMSE_MULTI_GLM_2 <- cbind(RMSE_MULTI_GLM_2,RMSE_MULTI_glm_2)
SD_MULTI_GLM_2 <- cbind(SD_MULTI_GLM_2,SD_MULTI_glm_2)


RMSE_MULTIP_GLM_2 <- RMSE_MULTI_GLM_2
SD_MULTIP_GLM_2 <- SD_MULTI_GLM_2 

RMSE_MULTIP_GLM_2 ; SD_MULTIP_GLM_2 