weightedAdaptive = function(y, K1,K2, X = NULL, r.all=seq(0, 1, 0.01)){
  n = length(y)
  if (is.null(X)) {
    p = 1
    X1 = matrix(1, nrow=n)
    mod = lm(y~1)#, family = "gaussian")
  } else {
    p = ncol(X)
    X1 = cbind(1, X)
    mod = lm(y~X)#, family = "gaussian")
  }
  s2 = summary(mod)$s**2  ## residual standard error ^2 , i.e sigma
  res = resid(mod)
  
  Z1.1<-K1 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K1)
  Z2.1<-K2 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K2)
  
  # Get individual p values 
  res.out = NULL
  
  ## Get Q
  out.Q<-SKAT_2Kernel_Optimal_Get_Q(K1, K2, res, r.all, n.Resampling = 0 , res.out)
  Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res)/s2
  
  ## Get p-values from Q
  out<-SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(Q.all, Z1.1/sqrt(2), 
                                             Z2.1/sqrt(2), r.all, method = "C")  
  
  param<-list(p.val.each=NULL,q.val.each=NULL)
  param$p.val.each<-out$p.val.each[1,]
  param$q.val.each<-Q.all[1,]
  param$rho<-r.all
  param$minp<-min(param$p.val.each)
  
  p.value<-out$p.value[1]
  
  re<-list(p.value = p.value, Test.Type = "C", Q = NA, param=param )  
  return(re)
  
}	

weightedLogistic = function(y, K1,K2, X = NULL, r.all=seq(0, 1, 0.01), pvalue){
  n = length(y)
  if (is.null(X)) {
    X1 <- rep(1, n)
    obj <- SKAT_Null_Model(y ~ 1, out_type = "D", Adjustment = F)
  } else {
    p = ncol(X)
    X1 = cbind(1, X)
    obj <- SKAT_Null_Model(y ~ X1, out_type = "D", Adjustment = F)
  }
  
  pi_1 <- obj$pi_1
  res <- obj$res
  
  Z1.1<- (K1 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (K1 * pi_1))
  Z2.1<- (K2 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (K2 * pi_1))
  
  s2<-1
  
  # Get individual p values 
  res.out = NULL
  
  ## Get Q
  out.Q<-SKAT_2Kernel_Optimal_Get_Q(K1, K2, res, r.all, n.Resampling = 0 , res.out)
  Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res)
  
  ## Get p-values from Q, something is up with the sqrt(2), change to 1.5
  out<-SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(Q.all, Z1.1/sqrt(2), 
                                             Z2.1/sqrt(2), r.all, method = "D")  
  
  if (pvalue == "Cauchy"){
    T <- mean(tan(pi*(0.5 - out$p.val.each)))
    p.val <- 0.5 - (atan(T))/pi
    
    return(p.val)
  }
  
  else{
    param<-list(p.val.each=NULL,q.val.each=NULL)
    param$p.val.each<-out$p.val.each[1,]
    param$q.val.each<-Q.all[1,]
    param$rho<-r.all
    param$minp<-min(param$p.val.each)
    
    p.value<-out$p.value[1]
    
    re<-list(p.value = p.value, Test.Type = "D", Q = NA, param=param )  
    return(re)
  }

}	

# For using WK
## First column of metabolite dataset must be KEGG
## metab.list contains pathway and metabolite connections
useWK <- function(y, metab, genomic, pathways, metab.list, X, type, pvalue = "NotCauchy"){
  n.pathways <- length(pathways)
  WK.pval <- matrix(rep(NA, n.pathways*2), nrow = n.pathways)
  for (i in c(1:n.pathways)){
    tryCatch({
      Z <- as.numeric(t(as.matrix(genomic[rownames(genomic) == pathways[i],])))
      W <- t(as.matrix(metab[rownames(metab)%in%metab.list[[pathways[i]]][,1],]))
      
      K_metab <- W%*%t(W)
      K.m_scale <- K_metab/sum(diag(K_metab), na.rm = T)
      K_genom <- Z%*%t(Z)
      K.g_scale <- K_genom/sum(diag(K_genom), na.rm = T)
      
      if (type == "D"){
        pval <- weightedLogistic(y, K.g_scale, K.m_scale, X = X, pvalue = pvalue)[[1]]
      }
      else{
        pval <- weightedAdaptive(y, K.g_scale, K.m_scale, X = X)[[1]]
      }
      
      WK.pval[i,] <- c(pathways[i], pval)},
      error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  WK.pval <- as.data.frame(WK.pval[complete.cases(WK.pval),])
  WK.pval$fdr <- p.adjust(as.numeric(as.character(WK.pval$V2)), method = "fdr")

  colnames(WK.pval) <- c("Pathway", "P-Value", "FDR")
  
  sig.WK <- WK.pval[WK.pval$FDR <= 0.05,]
  return(list(WK.pval, sig.WK))
}

useGenom <- function(y, genomic, pathways, metab.list, X, type){
  attach(X)
  n.pathways <- length(pathways)
  lin.pval <- matrix(rep(NA, n.pathways*2), nrow = n.pathways)
  
  for (i in c(1:length(pathways))){
    tryCatch({
      Z <- as.numeric(t(as.matrix(genomic[rownames(genomic) == pathways[i],])))
      
      var <- c("Z", colnames(X))
      model <- formula(paste("y", "~", paste(var, collapse=" + ")))
      if (type == "D"){
        mod <- glm(model, family = "binomial") 
      }
      else{
        mod <- lm(model)
      }
      
      lin.pval[i,] <- c(pathways[i], coef(summary(mod))[2,4])},
      error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  lin.pval <- as.data.frame(lin.pval[complete.cases(lin.pval),])
  lin.pval$fdr <- p.adjust(as.numeric(as.character(lin.pval$V2)), method = "fdr")
  
  colnames(lin.pval) <- c("Pathway", "P-Value", "FDR")

  sig.lin <- lin.pval[lin.pval$FDR <= 0.05,]
  return(list(lin.pval, sig.lin))
  
}

useMetab <- function(y, metab, pathways, type, X){
  n.pathways <- length(pathways)
  skat.pval <- matrix(rep(NA, n.pathways*2), nrow = n.pathways)
  
  if(is.null(X)){
    obj <- SKAT_Null_Model(y ~ 1, out_type = type, Adjustment = F)
  }
  else{
    obj <- SKAT_Null_Model(y ~ X, out_type = type, Adjustment = F) 
  }
  
  
  for (i in c(1:n.pathways)){
    tryCatch({
      metab1 <- t(as.matrix(metab[rownames(metab)%in%metab.list[[pathways[i]]][,1],]))
      pval <- SKAT(metab1, obj, kernel = "linear", is_check_genotype = F)$p.value
      skat.pval[i,] <- c(pathways[i], pval)},
      error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  skat.pval <- as.data.frame(skat.pval[complete.cases(skat.pval),])
  skat.pval$fdr <- p.adjust(as.numeric(as.character(skat.pval$V2)), method = "fdr")
  
  colnames(skat.pval) <- c("Pathway", "P-Value", "FDR")
  
  sig.skat <- skat.pval[skat.pval$FDR <= 0.05,]
  return(list(skat.pval, sig.skat))
}

weightedLongitudinal <- function(y, dat1, dat2, X = NULL, time, r.all=seq(0, 1, 0.01)){
  if (is.null(X)) {
    p = 1
    X1 = cbind(rep(1, n*t), id = rep(1:n, each = t))
    mod = lmer(y~1 + (1|X1[,2]))#, family = "gaussian")
  } else {
    p = ncol(X)
    time <- rep(1:3, n)
    X1 <- cbind(1, X)
    id = rep(1:n, each = t)
    #mod = lmer(y ~ X + time + (1|id), REML = T)
    mod = lmer(y~X + time + (time|id), REML = T)#, family = "gaussian") # Random slopes/intercept
  }
  res <- as.matrix(residuals(mod))
  
  V_est <- as.data.frame(VarCorr(mod)) # estimates for covariance matrix
  D_hat <- matrix(c(V_est[1, 4], V_est[3, 4], V_est[3, 4], V_est[2, 4]), nrow = 2)
  sigma_hat <- as.data.frame(VarCorr(mod))[4,4]
  y_var <- B%*%D_hat%*%t(B) # Random slope and intercept
  Sigma <- kronecker(diag(1, n), y_var) + diag(sigma_hat, n*t)
  sig_inv <- solve(Sigma)
  
  X1 <- cbind(1, X)
  
  K1 <- dat1%*%t(dat1)
  K1.scale <- K1/sum(diag(K1)) ## Molidify by the number of timepoints?
  
  K2 <- dat2%*%t(dat2)
  K2.scale <- K2/sum(diag(K2))
  
  Z1.1 <- K1.scale - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K1.scale)
  Z2.1 <- K2.scale - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K2.scale)
  
  # Get all Q
  Q.all <- rep(NA, length(r.all))
  
  for (i in 1:length(r.all)){
    K.rho <- r.all[i]*K1.scale + (1- r.all[i])*K2.scale
    Q.rho <- t(res)%*%sig_inv%*%K.rho%*%sig_inv%*%res
    Q.all[i] <- Q.rho
  }
  ## Try SKAT_2Kernel_Ortho_Optimal_Get_Pvalue from joint-common rare file!
  
  out<-SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(t(as.matrix(Q.all)), Z1.1/sqrt(2), 
                                             Z2.1/sqrt(2), r.all, method = "C")  
  
  
  ## Get all Q
  Q.all <- rep(NA, length(r.all))
  Q.pvalue <- rep(NA, length(r.all))
  P <- sig_inv - sig_inv%*%X1%*%solve(t(X1)%*%sig_inv%*%X1)%*%t(X1)%*%sig_inv
  p_decom <- eigen(P)
  phalf <- p_decom$vectors%*%diag(sqrt(round(p_decom$values, 6)))%*%t(p_decom$vectors)
  
  for (i in 1:length(r.all)){
    K.rho <- r.all[i]*K1.scale + (1- r.all[i])*K2.scale
    Q.rho <- t(res)%*%sig_inv%*%K.rho%*%sig_inv%*%res
    Q.all[i] <- Q.rho
    
    ## Get P square root through decomposition
    
    M <- r.all[i]*phalf%*%K1.scale%*%phalf + (1- r.all[i])*phalf%*%K2.scale%*%phalf
    eig <- eigen(M)
    values <- eig$values[eig$values>1e-6*eig$values[1]]
    
    Q.pvalue[i] <- davies(Q.rho, values)$Qq
    
    
  }
  
  Q.pvalue <- Q.pvalue[Q.pvalue < 1]
  T <- mean(tan(pi*(0.5 - Q.pvalue)))
  p.val <- 0.5 - (atan(T))/pi
  
  return(p.val)
  
}	

## For Cauchy combination test
weightedCauchy = function(y, K1,K2, X = NULL, r.all=seq(0, 1, 0.01)){
  n = length(y)
  if (is.null(X)) {
    p = 1
    X1 = matrix(1, nrow=n)
    mod = lm(y~1)#, family = "gaussian")
  } else {
    p = ncol(X)
    X1 = cbind(1, X)
    mod = lm(y~X)#, family = "gaussian")
  }
  s2 = summary(mod)$s**2  ## residual standard error ^2 , i.e sigma
  res = resid(mod)
  
  Z1.1<-K1 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K1)
  Z2.1<-K2 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K2)
  
  # Get individual p values 
  res.out = NULL
  
  ## Get Q's and all p-values from Q
  out.Q<-SKAT_2Kernel_Optimal_Get_Q(K1, K2, res, r.all, n.Resampling = 0 , res.out)
  Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res)/s2
  
  out<-SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(Q.all, Z1.1/sqrt(2), 
                                             Z2.1/sqrt(2), r.all, method = "C")
  

  ## Cauchy Combination
  
  T <- mean(tan(pi*(0.5 - out$p.val.each)))
  p.val <- 0.5 - (atan(T))/pi
  
  return(p.val)
}	

