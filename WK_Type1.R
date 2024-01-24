library(CompQuadForm)
library(mvtnorm)
args <- commandArgs(TRUE)
print(args)

load("RealData.rData")

n <- as.numeric(args[1])
pathway <- args[2]

# Commented out Kernel PCA part
try(source("Joint_CommonRare.R"))
try(source("WK_Scripts.R"))

doOne <- function(pathway){
  X <- cbind(rnorm(n, 0, 1), rbinom(n, 1, 0.5))
  X_tilde <- cbind(1, X)
  
  Z <- quantile(ecdf(as.numeric(pathway_values[pathway_values$value==pathway,][-1])), 
                runif(n, 0, 1))
  W.real <- t(metab.log[metab.log$Pathway==pathway, -c(1:2)])
  sigma <- cov(W.real)
  W <- rmvnorm(n, mean = rep(0, nrow(sigma)), sigma = sigma)
  
  #W.o <- (diag(1, n) - Z%*%solve(t(Z)%*%Z)%*%t(Z))%*%W
  alpha <- c(1, 0.5, 0.5)
  
  # Kernels, shouldn't need to scale due to weights
  K_metab <- W%*%t(W)
  K.m_scale <- K_metab/sum(diag(K_metab))
  K_genom <- Z%*%t(Z)
  K.g_scale <- K_genom/sum(diag(K_genom))
  
  epsilon <- rnorm(n, 0, 1)
  y <- X_tilde%*%alpha + epsilon 
  
  cauchy.p <- weightedCauchy(y, K.m_scale, K.g_scale, X)[[1]]
  min.p <- weightedAdaptive(y, K.m_scale, K.g_scale, X)[[1]]
  return(c(min.p, cauchy.p))
}

# Use real data, map02010 has 37 metabolites, map00564 has 6 metabolites

sim_05 <- t(replicate(5000, doOne(pathway)))
sim_01 <- t(replicate(5000, doOne(pathway)))
final <- data.frame(sim_05, sim_01)
colnames(final) <- c("a_05_minp", "a_05_cauchy", "a_01_minp", "a_01_cauchy")

filename <- paste("Type1_New/WK_n", n, pathway, ".csv", sep = "")
write.csv(final, filename)