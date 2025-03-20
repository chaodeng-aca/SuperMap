#' @title pairwise_optim
#' @description This function is used to solve the optimization problem corresponding to the pairwise interactions part in ADMM
#' @param Y The data matrix Y
#' @param X The data matrix X
#' @param B_initial The initial value of the coefficient matrix B
#' @param sigma_init The initial value of the variance sigma
#' @param index_B_01 The index matrix of non-zero positions in the coefficient matrix B
#' @param dist_matrix prior-related weight matrix D
#' @param lambda_B prior penalty parameter
#' @param rho dual variable parameter
#' @param w parameter for weighting marginal distributions and pairwise interactions
#' @param marginal_par The parameters obtained from solving the marginal distribution part in ADMM
#' @param u The dual variables updated in the last iteration
#' @param iteration The number of iterations
#' @param d0 The hyperparameter d0
#'
#' @return Return the solved parameters B, sigma, and the loss value for each iteration.
#' @export
#'
pairwise_optim = function(Y, X, B_initial, sigma_init,index_B_01, dist_matrix, lambda_B, rho, w, marginal_par, u,iteration,d0){
  
  R = (t(Y)%*%Y)/nrow(Y)
  r = (t(X)%*%X)/nrow(X)
  B = B_initial
  sigma_ele = sigma_init
  
  
  X_tilde_t = matrix(0,nrow = nrow(r),ncol = ncol(B))
  for (i in 1:ncol(X_tilde_t)) {
    if (sum(index_B_01[,i]==1)==1) {
      X_tilde_t[,i] = matrix(r[,index_B_01[,i]==1],ncol = 1) %*% B[index_B_01[,i]==1,i]
    }else{
      X_tilde_t[,i] = r[,index_B_01[,i]==1] %*% B[index_B_01[,i]==1,i]
    }
  }
  X_tilde<- t(X_tilde_t)
  
  if (iteration==1) {
    pairwise_iter = 6
  }else{
    pairwise_iter = 5
  }
  for (quadratic_iter in 1:pairwise_iter) {
    loss = pairwise_loss(R = R,r = r,B = B,sigma = sigma_ele,dist_matrix = dist_matrix/d0,marginal_par = marginal_par,u = u,lambda_prior = lambda_B,rho = rho,w = w)
    print(c(loss[[5 ]],loss[[4]],loss[[3]],loss[[1]],loss[[2]]))
    start_time = Sys.time()
    for (j in 1:ncol(B)) {
      index = (index_B_01[,j] == 1)
      
      prior = exp(dist_matrix[index,j]/d0)
      prior = c(prior,0)
      prior = diag(prior**2, nrow = length(prior),ncol = length(prior))
      
      e = rep(0,nrow(X_tilde))
      e[j] = 1
      
      if (sum(index)==1) {
        x_tmp = cbind(matrix(X_tilde[, index],ncol = 1),e)
        pairwise_res = pinv(t(x_tmp)%*%x_tmp + (rho/w)*diag(ncol(x_tmp))+ (lambda_B/w)*prior)%*%(t(x_tmp)%*%R[,j]+(rho/w)*(marginal_par[c(index,TRUE),j]+u[c(index,TRUE),j]))
        
        B[index,j] = pairwise_res[1:(nrow(pairwise_res)-1),]
        sigma_ele[j] = pairwise_res[nrow(pairwise_res),]
        X_tilde[j,]<- t(B[index,j]) %*% matrix(r[index,],nrow=1)
        
      }else{
        x_tmp = cbind(X_tilde[, index],e)
        pairwise_res = pinv(t(x_tmp)%*%x_tmp + (rho/w)*diag(ncol(x_tmp))+ (lambda_B/w)*prior)%*%(t(x_tmp)%*%R[,j]+(rho/w)*(marginal_par[c(index,TRUE),j]+u[c(index,TRUE),j]))
        B[index,j] = pairwise_res[1:(nrow(pairwise_res)-1),]
        sigma_ele[j] = pairwise_res[nrow(pairwise_res),]
        X_tilde[j,]<- t(B[index,j]) %*% r[index,]
      }
    }
    print(Sys.time()-start_time)
  }

  return(list(B,sigma_ele,loss[[1]],loss[[2]]))
}

pairwise_loss = function(R,r,B,sigma,dist_matrix,marginal_par,u,lambda_prior,rho,w){
  loss_pairwise = R - diag(sigma) - t(B)%*%r%*%B
  
  loss_penalty = (lambda_prior)*sum((exp(dist_matrix)*B)**2)
  loss_st = marginal_par - rbind(B,sigma) + u
  
  loss = w*sum(loss_pairwise^2) + loss_penalty + rho * sum(loss_st^2)
  loss = loss/2
  loss_g = w*sum(loss_pairwise^2)/2 + loss_penalty/2
  return(list(loss,loss_g,loss_penalty/2,rho * sum(loss_st^2)/2,sum((marginal_par - rbind(B,sigma))^2)))
}


#' @title nls_diag
#' @description This function utilizes the Levenberg-Marquardt (LM) method to solve the nonlinear least squares problem
#' @export
nls_diag = function(Y, X, pairwise_par, marginal_init, u,rho){
  m<- length(Y)
  Y_unpaired_sort<- sort(Y)    
  density_Ym = (1:m)/m
  
  
  f_xy_function<- function(par,Y_sort_ele,density_Y_ele){
    f_x<- function(y){
      f_x<- 0
      for(i in 1:nrow(X)){
        f_x<- f_x + 1/nrow(X) * 1/(1+exp((X[i,]%*%par[-length(par)]-y)*pi/(3*par[length(par)])^0.5))
      }
      return(f_x)
    }
    f_xy<- ((density_Y_ele - f_x(Y_sort_ele))^2 + (rho/2)/m * sum((par - pairwise_par + u)^2))^0.5
    return(f_xy)
  }

  lower_bound = rep(-Inf,length(pairwise_par))
  lower_bound[length(lower_bound)] = 10^(-10)  
  upper_bound = rep(Inf,length(pairwise_par))
  upper_bound[length(upper_bound)] = var(Y)
  
  model<- nls.lm(par = marginal_init, lower = lower_bound, upper = upper_bound, fn = f_xy_function, Y_sort_ele = Y_unpaired_sort, density_Y_ele = density_Ym,
                 control = nls.lm.control(maxiter=300))
  return(model)
}



loss_f = function(Y, X, parameter){
  m<- length(Y)
  Y_unpaired_sort<- sort(Y)     
  density_Ym<- (1:m)/m
  f_x<- function(y){
    f_x<- 0
    for(i in 1:nrow(X)){
      #f_x<- f_x + 1/nrow(X) * 1/(1+exp((X[i,]%*%parameter[-length(parameter)]-y)/parameter[length(parameter)]))
      f_x<- f_x + 1/nrow(X) * 1/(1+exp((X[i,]%*%parameter[-length(parameter)]-y)*pi/(3*parameter[length(parameter)])^0.5))
    }
    return(f_x)
  }
  f_xy<- sum((density_Ym - f_x(Y_unpaired_sort))^2)
  return(f_xy)
  
}

#' @title marginal_nls
#' @description This function is used to solve the optimization problem corresponding to the marginal distribution part in ADMM
#' @param Y The data matrix Y
#' @param X The data matrix X
#' @param pairwise_par The parameters obtained from solving the pairwise interactions part in ADMM
#' @param marginal_init_des The initial values of the parameters to be solved
#' @param u The dual variables updated in the last iteration
#' @param index_B_01_des The index matrix of non-zero positions in the coefficient matrix B
#' @param nls_diag The custom non-negative least squares solver function that needs to be called
#' @param loss_f The custom loss function
#' @param rho dual variable parameter
#'
#' @return The solved parameters and the loss function values
#' @export

marginal_nls = function(Y, X, pairwise_par, marginal_init_des, u, index_B_01_des, nls_diag, loss_f, rho,n){
  foreach(j = 1:n, .combine = rbind, .errorhandling = "pass",.packages = c('bigmemory','minpack.lm'),.export=c('Y', 'X', 'pairwise_par', 'u','index_B_01_des')) %dopar% {
    Y_unpaired <- attach.big.matrix(Y)
    X_unpaired <- attach.big.matrix(X)
    pairwise_parameter <- attach.big.matrix(pairwise_par)
    marginal_init = attach.big.matrix(marginal_init_des)
    u_iter = attach.big.matrix(u)
    index_B_01 <- attach.big.matrix(index_B_01_des)
    
    if (sum(index_B_01[,j]==1)==1) {  
      model_nls <- nls_diag(Y = Y_unpaired[,j],X = matrix(X_unpaired[,index_B_01[,j] == 1],ncol = 1),
                            pairwise_par = pairwise_parameter[c(index_B_01[,j]==1,TRUE),j], marginal_init = marginal_init[c(index_B_01[,j]==1,TRUE),j],
                            u = u_iter[c(index_B_01[,j]==1,TRUE),j],rho = rho)
      loss_value_f = loss_f(Y = Y_unpaired[,j],X = matrix(X_unpaired[,index_B_01[,j] == 1],ncol = 1),parameter = model_nls$par)
    }else{
      model_nls <- nls_diag(Y = Y_unpaired[,j],X = X_unpaired[,index_B_01[,j] == 1],
                            pairwise_par = pairwise_parameter[c(index_B_01[,j]==1,TRUE),j], marginal_init = marginal_init[c(index_B_01[,j]==1,TRUE),j],
                            u = u_iter[c(index_B_01[,j]==1,TRUE),j],rho = rho)
      loss_value_f = loss_f(Y = Y_unpaired[,j],X = X_unpaired[,index_B_01[,j] == 1],parameter = model_nls$par)
    }
    
    loss_diag = model_nls$deviance
    convergence = model_nls$info
    marginal_par = model_nls$par
    
    marginal_alpha = rep(0,nrow(index_B_01)+1)
    marginal_alpha[c(index_B_01[,j] == 1,TRUE)] = marginal_par
    
  
    c(marginal_alpha,loss_diag,convergence,loss_value_f)
  }
}



#' @title supermap
#' @description This function utilizes the ADMM algorithm to estimate parameters in SuperMap model.
#' @param Y_unpaired The data matrix Y
#' @param X_unpaired The data matrix X
#' @param index_B_nonzero The index matrix of non-zero positions in the coefficient matrix B
#' @param B_dist prior-related weight matrix D
#' @param w parameter for weighting marginal distributions and pairwise interactions
#' @param rho dual variable parameter
#' @param lambda_prior prior penalty parameter
#' @param n_iter The number of iterations
#' @param d0 The hyperparameter d0, with a default value of 100,000
#'
#' @return The estimated model parameters, as well as the loss value for each iteration
#' @export
#'
supermap = function(Y_unpaired, X_unpaired, index_B_nonzero, B_dist, w, rho, lambda_prior, n_iter=20,d0=100000,ncore=50){
  
  
  X_mean = colMeans(X_unpaired)
  Y_mean = colMeans(Y_unpaired)
  X_unpaired = apply(X_unpaired, 2, function(x) x-mean(x))
  Y_unpaired = apply(Y_unpaired, 2, function(x) x-mean(x))
  
  
  bigmat_Y_unpaired = as.big.matrix(Y_unpaired)
  bigmat_X_unpaired = as.big.matrix(X_unpaired)
  bigmat_index_B_01 = as.big.matrix(index_B_nonzero)
  descriptor_Y <- describe(bigmat_Y_unpaired)
  descriptor_X <- describe(bigmat_X_unpaired)
  descriptor_index <- describe(bigmat_index_B_01)
  
  Y_col_variances = apply(Y_unpaired,2,var)
  B_pairwise = index_B_nonzero*0.5
  sigma_pairwise = 0.5*Y_col_variances
  marginal_par = rbind(index_B_nonzero*0.5,0.5*Y_col_variances)
  u_init = matrix(0,nrow = nrow(marginal_par),ncol = ncol(marginal_par))
  
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  
  convergence_matrix = matrix(0,nrow = n_iter,ncol = 2)
  gene_cor_matrix = matrix(nrow = n_iter,ncol = ncol(Y_unpaired))
  cov_cor_matrix = matrix(nrow = n_iter,ncol = ncol(Y_unpaired))
  
  for (iteration in 1:n_iter) {
    print(iteration)
    pairwise_result = pairwise_optim(Y = Y_unpaired,X = X_unpaired,B_initial = B_pairwise,sigma_init = sigma_pairwise,index_B_01 = index_B_nonzero,dist_matrix = B_dist,
                                     lambda_B = lambda_prior,rho = rho,w = w,marginal_par = marginal_par,u = u_init,iteration = iteration,d0 = d0)
    pairwise_par = rbind(pairwise_result[[1]],pairwise_result[[2]])
    
    bigmat_pairwise_par = as.big.matrix(pairwise_par)
    bigmat_marginal_par = as.big.matrix(marginal_par)
    bigmat_u = as.big.matrix(u_init)
    descriptor_pairwise_par <- describe(bigmat_pairwise_par)
    descriptor_marginal_par <- describe(bigmat_marginal_par)
    descriptor_u <- describe(bigmat_u)
    start_time = Sys.time()
    # marginal_result = marginal_optim(Y = descriptor_Y, X = descriptor_X, pairwise_par = descriptor_pairwise_par,marginal_init_des = descriptor_marginal_par,
    #                                  u = descriptor_u,index_B_01_des = descriptor_index, optim_diag = optim_diag,n=ncol(Y_unpaired),rho = rho,loss_f = loss_f)
    marginal_result = marginal_nls(Y = descriptor_Y, X = descriptor_X, pairwise_par = descriptor_pairwise_par,marginal_init_des = descriptor_marginal_par,
                                   u = descriptor_u,index_B_01_des = descriptor_index, nls_diag = nls_diag,n=ncol(Y_unpaired),rho = rho,loss_f = loss_f)
    print(Sys.time()-start_time)
    marginal_result = t(marginal_result)
    marginal_par = marginal_result[1:nrow(pairwise_par),]
    marginal_loss = marginal_result[nrow(pairwise_par)+1,]
    marginal_convergence = marginal_result[nrow(pairwise_par)+2,]
    
    
    u_init = u_init + marginal_par - pairwise_par
    B_pairwise = pairwise_result[[1]]
    sigma_pairwise = pairwise_result[[2]]
    
    r = sum((marginal_par - pairwise_par)^2)
    g_loss = pairwise_result[[4]]
    f_loss = sum(marginal_result[nrow(pairwise_par)+3,])
    cat(g_loss,f_loss,g_loss+f_loss,r,'\n')
    convergence_matrix[iteration,1] = r
    convergence_matrix[iteration,2] = g_loss+f_loss
    
    
  }
  intercept_estimated = as.vector(Y_mean - t(X_mean)%*%B_pairwise)
  
  stopImplicitCluster()
  stopCluster(cl)
  return(list(estimate_b = B_pairwise,estimated_intercept = intercept_estimated,estimated_sigma = pairwise_result[[2]],
              convergence=convergence_matrix))
}



#' @title getmetacell_RNA
#' @description This function allows us to obtain the gene expression matrix of metacells
#' @param object A SeuratObject 
#' @param cluster.name The cluster name in the SeuratObject
#' @param size size factor for normalization
#'
#' @return A gene expression matrix, where rows represent cells and columns represent genes
#' @export
getmetacell_RNA <- function(object, cluster.name, size_factor=10^6) {
  clusters <- object@meta.data[[cluster.name]]
  clust.levels <- levels(clusters)
  assay.matrix = object@assays$RNA$counts
  nrow <- length(clust.levels)
  ncol <- nrow(assay.matrix)
  metacell <- matrix(nrow = nrow, ncol = ncol)
  rownames(metacell) <- paste0("metacell_", clust.levels)
  colnames(metacell) <- rownames(assay.matrix)
  for (i in 1:nrow(metacell)) {
    metacell_expression <- apply(assay.matrix[ ,clusters ==as.character(i - 1)], 1, mean)
    metacell[i,] = log(metacell_expression*size_factor/sum(metacell_expression)+1)
  }
  return(metacell)
}

#' @title getmetacell_ATAC
#' @description This function allows us to obtain the chromatin accessibility matrix of metacells
#' 
#' @param object A SeuratObject
#' @param cluster.name The cluster name in the SeuratObject
#' @return  A chromatin accessibility matrix, where rows represent cells and columns represent peaks
#' @export
getmetacell_ATAC <- function(object, cluster.name) {
  clusters <- object@meta.data[[cluster.name]]
  clust.levels <- levels(clusters)
  assay.matrix = object@assays$ATAC@counts
  assay.matrix = ifelse(as.matrix(assay.matrix) == 0, 0, 1)
  assay.matrix = Matrix(assay.matrix,sparse = TRUE)
  nrow <- length(clust.levels)
  ncol <- nrow(assay.matrix)
  metacell <- matrix(nrow = nrow, ncol = ncol)
  rownames(metacell) <- paste0("metacell_", clust.levels)
  colnames(metacell) <- rownames(assay.matrix)
  for (i in 1:nrow(metacell)) {
    metacell[i, ] <- apply(assay.matrix[, clusters ==as.character(i - 1)], 1, mean)
  }
  return(metacell)
}



