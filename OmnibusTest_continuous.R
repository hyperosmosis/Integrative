nperm = 1000
YYY = matrix(rnorm(nperm*2500), ncol = nperm)  # assume that there are no more than 100 SNPs for this simulation 
# Only do 1000 permutations 
OmnibusTest_Continuous = function(y, Ks, X = NULL, nperm = 1000) {
  try(library(CompQuadForm))
  try(library(CompQuadForm, lib.loc = "~/Rlib"))
  
  ## Helper Functions
  getDavies = function(Qs, evs) {
    ## function for getting the davies p-values
    return(davies(Qs, evs)$Qq)
  }
  
  getIndivP = function(K, M, A, s2) {
    ## gets all of the parameters necessary for each kernel
    S = t(M)%*%K%*%M/s2
    W = A%*%K%*%A
    scale = sum(diag(W))
    ee = eigen(W, symmetric = T)
    w = which(ee$values>1e-10)
    pvals = davies(S,ee$values[w])$Qq
    return(list(S=S, W = W,scale = scale, m = length(w), evals = ee$values[w], evecs = ee$vectors[,w], pvals = pvals))
  }
  
  ## This is just temporary stuff for testing code
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
  s2 = summary(mod)$s**2      ## residual standard error ^2 , i.e sigma
  
  D0 = diag(n)#diag(mu*(1-mu))
  M = resid(mod)
  ##P0= D0 - D0%*% X1%*%solve(t(X1)%*%D0%*%X1)%*%t(X1)%*%D0
  P0= D0 - X1%*%solve(t(X1)%*%X1)%*%t(X1)
  
  A = P0# V%*%sqrt(D)%*%t(V)
  ## Get individual p-values
  S = sapply(Ks, getIndivP, M, A, s2)
  pvals = unlist(S[7,])
  Q.stat = min(pvals)
  ms = unlist(S[4,])
  V_vector = unlist(S[6,])
  V = matrix(V_vector, n, sum(ms), byrow = F)
  Sigma = t(V) %*% V
  
  # This implies that the first block will be identity matrix
  
  e = eigen(t(Sigma), symmetric = T)
  w = which(e$values>1e-10)
  # Is it possible to have that w to be zero?
  e$values[-w] = 0
  L =  t(sqrt(e$values)*t(e$vectors))
  rstar = (L%*%YYY[1:sum(ms),])**2
  
  plat = matrix(nrow = length(Ks), ncol = nperm)
  
  for (i in seq(length(Ks))){
    if (i == 1) {
      rows = seq(ms[1])
    } else {
      rows = seq(ms[i])+sum(ms[seq(i-1)])
    }
    if (length(rows) > 1){
      plat[i,] = sapply(colSums(rstar[rows,]*S[,i]$evals),getDavies, S[,i]$evals)
    }else if (length(rows) == 1){
      plat[i,] = sapply(rstar[rows,]*S[,i]$evals,getDavies, S[,i]$evals)
    }
  }
  Q.perm = apply(plat,2, min) # These were actually p values in this case
  resamp.pval = mean(Q.perm<=Q.stat)
  return(list(indiv.pvals = pvals, omnibus.pval = resamp.pval))
} 

# library(SKAT, lib ="/netscr/nzhao/Dissertation/Rlib")   # This was the old version of SKAT 
# Obtained function from SKAT.0.91 
# try(source("/netscr/nzhao/Dissertation/Simulation/OmnibusTest/Joint_CommonRare.R"))
try(source("Joint_CommonRare.R"))
OmnibusTest_adaptive = function(y, K1,K2, X = NULL, r.all=c(0,0.25,0.5,0.75,1)){ 
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
  s2 = summary(mod)$s**2      ## residual standard error ^2 , i.e sigma
  #D0 = diag(n)#diag(mu*(1-mu))
  res = resid(mod)
  # P0= D0 - X1%*%solve(t(X1)%*%X1)%*%t(X1)

  # Kernel PCA  to get the basis
  egK1 = eigen(K1, symmetric = T)
  #Z1 = egK1$vectors[,egK1$values > 1e-10] %*% diag(sqrt(egK1$values[egK1$values > 1e-10]))
  egK2 = eigen(K2, symmetric = T)
  #Z2 = egK2$vectors[,egK2$values > 1e-10] %*% diag(sqrt(egK2$values[egK2$values > 1e-10]))
  # project Z2 onto Z1 and obtain perpendicular part
  #out.QR<-qr(Z1)

	#Z1.Q<-try(as.matrix(qr.Q(out.QR)), silent = TRUE)  
	#Z1.R<-try(as.matrix(qr.R(out.QR)), silent = TRUE)
	
  #Z2.item<-t(Z1.Q) %*% Z2
	#Z2.item1<- Z1.Q %*% Z2.item
	#Z2.Ortho<-Z2 - Z2.item1
	#Z2<-Z2.Ortho

  Z1.1<-K1 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K1)
	Z2.1<-K2 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% K2)
  
  # Compute the Q values 
 	n.r<-length(r.all)
	p.m<-dim(K1)[2]
	Q.r<-rep(0,n.r)
  
  temp1<-t(mod$res) %*% K1
	temp2<-t(mod$res) %*% K2
	for(i in 1:n.r){
		r.corr<-r.all[i]
		Q1<-(1-r.corr) * rowSums(temp1^2)
		Q2<-r.corr * rowSums(temp2^2)
		Q.r[i]<-Q1 + Q2
	}
	Q.r = Q.r/(2*s2)   # Keep track of this especially 
 
  Q.all =list(Q.r=Q.r, Q.r.res=NULL , Q.sim=NULL)
	# Get individual p values 
  res.out = NULL
  
	out.Q<-SKAT_2Kernel_Optimal_Get_Q(K1, K2, res, r.all, n.Resampling = 0 , res.out)
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) / s2
	
  out<-SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(Q.all, Z1.1/sqrt(2), 
    Z2.1/sqrt(2), r.all, method = "C")  
    
 	param<-list(p.val.each=NULL,q.val.each=NULL)
	param$p.val.each<-out$p.val.each[1,]
	param$q.val.each<-Q.all[1,]
	param$rho<-r.all
	param$minp<-min(param$p.val.each)

	id_temp<-which(param$p.val.each == min(param$p.val.each))
	id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
	if(length(id_temp1) > 0){
		param$rho[id_temp1] = 1
	}
	param$rho_est<-param$rho[id_temp]


	p.value<-out$p.value[1]
  n.Resampling = 0 
	p.value.resampling<-NULL
	if(n.Resampling > 0){
		p.value.resampling<-out$p.value[-1]
		param$pval.each.resample<-out$p.val.each[-1,]
	}

 	re<-list(p.value = p.value, Test.Type = "C", Q = NA, param=param )  
  return(re)
  }	
  
  # May consider resampling approach later on 

