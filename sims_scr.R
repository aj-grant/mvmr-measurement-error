################################################################################
#Script to re-produce the primary simulation study
################################################################################

library(MendelianRandomization)
library(parallel)
library(MVMR)
library(tidyverse)
library(RColorBrewer)

N = 20000
p = 40
q = qnorm(0.975, 0, 1)

M = 1000

cl = makeCluster(6)
clusterEvalQ(cl, library(MendelianRandomization))
clusterEvalQ(cl, library(MVMR))
clusterExport(cl, c('N', 'p', 'M', 'sstat', 'fcalc'))

theta1 = 0.2
theta2 = 0
rho = c(0, 0.2, 0.4, 0.6, 0.8)
err = seq(0, 4, length = 9)
clusterExport(cl, c('theta1', 'theta2'))

#X1 causal, X2 null
clusterSetRNGStream(cl, 20210525)
D1 = lapply(1:5, function(i){
  lapply(1:9, function(j){
    rr = rho[i]
    ee = sqrt(err[j])
    parLapply(cl, 1:M, function(k){
      g = 0.2
      bx1 = runif(p, 0.08, 0.2)
      bx2 = runif(p, 0.08, 0.2)
      U = rnorm(2 * N, 0, 1)
      maf = runif(p, 0.01, 0.5)
      G = sapply(1:p, function(l){rbinom(2 * N, 2, maf[l])})
      e1 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
      e2 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
      x1 = G %*% bx1 + e1
      x2 = rr * x1 + sqrt(1 - rr^2) * (G %*% bx2 + e2)
      zeta1 = ee * rnorm(2 * N, 0, 1)
      X1 = x1 + zeta1
      X2 = x2
      Y = theta1 * x1 + theta2 * x2 + U + rnorm(2 * N, 0, 1)
      bxhat1 = bxhat2 = byhat = vector(length = p)
      sebx1 = sebx2 = seby = vector(length = p)
      sam1 = 1:N
      sam2 = (N+1):(2*N)
      for (l in 1:p){
        xmod1 = sstat(X1[sam1], G[sam1, l])
        bxhat1[l] = xmod1$bhat
        sebx1[l] = xmod1$se
        xmod2 = sstat(X2[sam1], G[sam1, l])
        bxhat2[l] = xmod2$bhat
        sebx2[l] = xmod2$se
        ymod = sstat(Y[sam2], G[sam2, l])
        byhat[l] = ymod$bhat
        seby[l] = ymod$se
      }
      fstat = fcalc(X1[sam1], X2[sam1], G[sam1, ])
      corX = cor(cbind(X1, X2))
      return(list("bxhat1" = bxhat1, "bxhat2" = bxhat2, "sebx1" = sebx1, "sebx2" = sebx2,
                  "byhat" = byhat, "seby" = seby,
                  "corbx1bx2" = cor(bxhat1, bxhat2), "varbx1" = var(bxhat1), "varbx2" = var(bxhat2),
                  "fstat" = fstat, "corX" = corX))
    })
  })
})

#X2 causal, X1 null
theta1 = 0
theta2 = 0.2
clusterExport(cl, c('theta1', 'theta2'))
clusterSetRNGStream(cl, 20210525)
D2 = lapply(1:5, function(i){
  lapply(1:9, function(j){
    rr = rho[i]
    ee = sqrt(err[j])
    parLapply(cl, 1:M, function(k){
      g = 0.2
      bx1 = runif(p, 0.08, 0.2)
      bx2 = runif(p, 0.08, 0.2)
      U = rnorm(2 * N, 0, 1)
      maf = runif(p, 0.01, 0.5)
      G = sapply(1:p, function(l){rbinom(2 * N, 2, maf[l])})
      e1 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
      e2 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
      x1 = G %*% bx1 + e1
      x2 = rr * x1 + sqrt(1 - rr^2) * (G %*% bx2 + e2)
      zeta1 = ee * rnorm(2 * N, 0, 1)
      X1 = x1 + zeta1
      X2 = x2
      Y = theta1 * x1 + theta2 * x2 + U + rnorm(2 * N, 0, 1)
      bxhat1 = bxhat2 = byhat = vector(length = p)
      sebx1 = sebx2 = seby = vector(length = p)
      sam1 = 1:N
      sam2 = (N+1):(2*N)
      for (l in 1:p){
        xmod1 = sstat(X1[sam1], G[sam1, l])
        bxhat1[l] = xmod1$bhat
        sebx1[l] = xmod1$se
        xmod2 = sstat(X2[sam1], G[sam1, l])
        bxhat2[l] = xmod2$bhat
        sebx2[l] = xmod2$se
        ymod = sstat(Y[sam2], G[sam2, l])
        byhat[l] = ymod$bhat
        seby[l] = ymod$se
      }
      fstat = fcalc(X1[sam1], X2[sam1], G[sam1, ])
      corX = cor(cbind(X1, X2))
      return(list("bxhat1" = bxhat1, "bxhat2" = bxhat2, "sebx1" = sebx1, "sebx2" = sebx2,
                  "byhat" = byhat, "seby" = seby,
                  "corbx1bx2" = cor(bxhat1, bxhat2), "varbx1" = var(bxhat1), "varbx2" = var(bxhat2),
                  "fstat" = fstat, "corX" = corX))
    })
  })
})

#X1 and X2 causal, X2 measured with error
theta1 = 0.2
theta2 = 0.2
clusterExport(cl, c('theta1', 'theta2'))
clusterSetRNGStream(cl, 20210525)
D3 = lapply(1:5, function(i){
  lapply(1:9, function(j){
    rr = rho[i]
    ee = sqrt(err[j])
    parLapply(cl, 1:M, function(k){
      g = 0.2
      bx1 = runif(p, 0.08, 0.2)
      bx2 = runif(p, 0.08, 0.2)
      U = rnorm(2 * N, 0, 1)
      maf = runif(p, 0.01, 0.5)
      G = sapply(1:p, function(l){rbinom(2 * N, 2, maf[l])})
      e1 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
      e2 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
      x1 = G %*% bx1 + e1
      x2 = rr * x1 + sqrt(1 - rr^2) * (G %*% bx2 + e2)
      zeta1 = ee * rnorm(2 * N, 0, 1)
      zeta2 = rnorm(2 * N, 0, 1)
      X1 = x1 + zeta1
      X2 = x2 + zeta2
      Y = theta1 * x1 + theta2 * x2 + U + rnorm(2 * N, 0, 1)
      bxhat1 = bxhat2 = byhat = vector(length = p)
      sebx1 = sebx2 = seby = vector(length = p)
      sam1 = 1:N
      sam2 = (N+1):(2*N)
      for (l in 1:p){
        xmod1 = sstat(X1[sam1], G[sam1, l])
        bxhat1[l] = xmod1$bhat
        sebx1[l] = xmod1$se
        xmod2 = sstat(X2[sam1], G[sam1, l])
        bxhat2[l] = xmod2$bhat
        sebx2[l] = xmod2$se
        ymod = sstat(Y[sam2], G[sam2, l])
        byhat[l] = ymod$bhat
        seby[l] = ymod$se
      }
      fstat = fcalc(X1[sam1], X2[sam1], G[sam1, ])
      corX = cor(cbind(X1, X2))
      return(list("bxhat1" = bxhat1, "bxhat2" = bxhat2, "sebx1" = sebx1, "sebx2" = sebx2,
                  "byhat" = byhat, "seby" = seby,
                  "corbx1bx2" = cor(bxhat1, bxhat2), "varbx1" = var(bxhat1), "varbx2" = var(bxhat2),
                  "fstat" = fstat, "corX" = corX))
    })
  })
})
stopCluster(cl)

set.seed(20210601)
#Estimates on D1
theta1 = 0.2
theta2 = 0
R1 = lapply(1:5, function(i){
  lapply(1:9, function(j){
    B = sapply(1:M, function(k){
      D = D1[[i]][[j]][[k]]
      mrob = mr_mvinput(bx = cbind(D$bxhat1, D$bxhat2), bxse = cbind(D$sebx1, D$sebx2), by = D$byhat, byse = D$seby)
      mrest = mr_mvivw(mrob)
      mrest_me = mrest_me(mrob, no_ini = 1)
      mrest_me_cor = mrest_me_cor(mrob, D$corX, no_ini = 1)
      c(mrest$Estimate, mrest$StdError, mrest_me$thest, sqrt(mrest_me$Var[1, 1]), sqrt(mrest_me$Var[2, 2]),
        mrest_me_cor$thest, sqrt(mrest_me_cor$Var[1, 1]), sqrt(mrest_me_cor$Var[2, 2]), D$fstat)
    })
  })
})

#Estimates on D2
theta1 = 0
theta2 = 0.2
R2 = lapply(1:5, function(i){
  lapply(1:9, function(j){
    B = sapply(1:M, function(k){
      D = D2[[i]][[j]][[k]]
      mrob = mr_mvinput(bx = cbind(D$bxhat1, D$bxhat2), bxse = cbind(D$sebx1, D$sebx2), by = D$byhat, byse = D$seby)
      mrest = mr_mvivw(mrob)
      mrest_me = mrest_me(mrob, no_ini = 1)
      mrest_me_cor = mrest_me_cor(mrob, D$corX, no_ini = 1)
      c(mrest$Estimate, mrest$StdError, mrest_me$thest, sqrt(mrest_me$Var[1, 1]), sqrt(mrest_me$Var[2, 2]),
        mrest_me_cor$thest, sqrt(mrest_me_cor$Var[1, 1]), sqrt(mrest_me_cor$Var[2, 2]), D$fstat)
    })
  })
})

#Estimates on D3
theta1 = 0.2
theta2 = 0.2
R3 = lapply(1:5, function(i){
  lapply(1:9, function(j){
    B = sapply(1:M, function(k){
      D = D3[[i]][[j]][[k]]
      mrob = mr_mvinput(bx = cbind(D$bxhat1, D$bxhat2), bxse = cbind(D$sebx1, D$sebx2), by = D$byhat, byse = D$seby)
      mrest = mr_mvivw(mrob)
      mrest_me = mrest_me(mrob, no_ini = 1)
      mrest_me_cor = mrest_me_cor(mrob, D$corX, no_ini = 1)
      c(mrest$Estimate, mrest$StdError, mrest_me$thest, sqrt(mrest_me$Var[1, 1]), sqrt(mrest_me$Var[2, 2]),
        mrest_me_cor$thest, sqrt(mrest_me_cor$Var[1, 1]), sqrt(mrest_me_cor$Var[2, 2]), D$fstat)
    })
  })
})

########################################################################################################################
#Tables
########################################################################################################################
theta1 = 0.2
theta2 = 0
T1 = matrix(nrow = 45, ncol = 36)
for (i in 1:5){
  for (j in 1:9){
    T1[(9*(i-1)+j), ] = c(rho[i], err[j],
                          mean(R1[[i]][[j]][1, ]), sd(R1[[i]][[j]][1, ]), mean(R1[[i]][[j]][3, ]), sum(R1[[i]][[j]][1, ] - q * R1[[i]][[j]][3, ] < theta1 & R1[[i]][[j]][1, ] + q * R1[[i]][[j]][3, ] > theta1)/M, sum(abs(R1[[i]][[j]][1, ]) - q * R1[[i]][[j]][3, ] > 0)/M,
                          mean(R1[[i]][[j]][2, ]), sd(R1[[i]][[j]][2, ]), mean(R1[[i]][[j]][4, ]), sum(R1[[i]][[j]][2, ] - q * R1[[i]][[j]][4, ] < theta2 & R1[[i]][[j]][2, ] + q * R1[[i]][[j]][4, ] > theta2)/M, sum(abs(R1[[i]][[j]][2, ]) - q * R1[[i]][[j]][4, ] > 0)/M,
                          
                          mean(R1[[i]][[j]][5, ]), sd(R1[[i]][[j]][5, ]), mean(R1[[i]][[j]][7, ]), sum(R1[[i]][[j]][5, ] - q * R1[[i]][[j]][7, ] < theta1 & R1[[i]][[j]][5, ] + q * R1[[i]][[j]][7, ] > theta1)/M, sum(abs(R1[[i]][[j]][5, ]) - q * R1[[i]][[j]][7, ] > 0)/M,
                          mean(R1[[i]][[j]][6, ]), sd(R1[[i]][[j]][6, ]), mean(R1[[i]][[j]][8, ]), sum(R1[[i]][[j]][6, ] - q * R1[[i]][[j]][8, ] < theta2 & R1[[i]][[j]][6, ] + q * R1[[i]][[j]][8, ] > theta2)/M, sum(abs(R1[[i]][[j]][6, ]) - q * R1[[i]][[j]][8, ] > 0)/M,
                          
                          mean(R1[[i]][[j]][9, ]), sd(R1[[i]][[j]][9, ]), mean(R1[[i]][[j]][11, ]), sum(R1[[i]][[j]][8, ] - q * R1[[i]][[j]][11, ] < theta1 & R1[[i]][[j]][9, ] + q * R1[[i]][[j]][11, ] > theta1)/M, sum(abs(R1[[i]][[j]][9, ]) - q * R1[[i]][[j]][11, ] > 0)/M,
                          mean(R1[[i]][[j]][10, ]), sd(R1[[i]][[j]][10, ]), mean(R1[[i]][[j]][12, ]), sum(R1[[i]][[j]][10, ] - q * R1[[i]][[j]][12, ] < theta2 & R1[[i]][[j]][10, ] + q * R1[[i]][[j]][12, ] > theta2)/M, sum(abs(R1[[i]][[j]][10, ]) - q * R1[[i]][[j]][12, ] > 0)/M,
                          
                          rowMeans(R1[[i]][[j]][13:16, ]))
  }
}

theta1 = 0
theta2 = 0.2
T2 = matrix(nrow = 45, ncol = 36)
for (i in 1:5){
  for (j in 1:9){
    T2[(9*(i-1)+j), ] = c(rho[i], err[j],
                          mean(R2[[i]][[j]][1, ]), sd(R2[[i]][[j]][1, ]), mean(R2[[i]][[j]][3, ]), sum(R2[[i]][[j]][1, ] - q * R2[[i]][[j]][3, ] < theta1 & R2[[i]][[j]][1, ] + q * R2[[i]][[j]][3, ] > theta1)/M, sum(abs(R2[[i]][[j]][1, ]) - q * R2[[i]][[j]][3, ] > 0)/M,
                          mean(R2[[i]][[j]][2, ]), sd(R2[[i]][[j]][2, ]), mean(R2[[i]][[j]][4, ]), sum(R2[[i]][[j]][2, ] - q * R2[[i]][[j]][4, ] < theta2 & R2[[i]][[j]][2, ] + q * R2[[i]][[j]][4, ] > theta2)/M, sum(abs(R2[[i]][[j]][2, ]) - q * R2[[i]][[j]][4, ] > 0)/M,
                          
                          mean(R2[[i]][[j]][5, ]), sd(R2[[i]][[j]][5, ]), mean(R2[[i]][[j]][7, ]), sum(R2[[i]][[j]][5, ] - q * R2[[i]][[j]][7, ] < theta1 & R2[[i]][[j]][5, ] + q * R2[[i]][[j]][7, ] > theta1)/M, sum(abs(R2[[i]][[j]][5, ]) - q * R2[[i]][[j]][7, ] > 0)/M,
                          mean(R2[[i]][[j]][6, ]), sd(R2[[i]][[j]][6, ]), mean(R2[[i]][[j]][8, ]), sum(R2[[i]][[j]][6, ] - q * R2[[i]][[j]][8, ] < theta2 & R2[[i]][[j]][6, ] + q * R2[[i]][[j]][8, ] > theta2)/M, sum(abs(R2[[i]][[j]][6, ]) - q * R2[[i]][[j]][8, ] > 0)/M,
                          
                          mean(R2[[i]][[j]][9, ]), sd(R2[[i]][[j]][9, ]), mean(R2[[i]][[j]][11, ]), sum(R2[[i]][[j]][8, ] - q * R2[[i]][[j]][11, ] < theta1 & R2[[i]][[j]][9, ] + q * R2[[i]][[j]][11, ] > theta1)/M, sum(abs(R2[[i]][[j]][9, ]) - q * R2[[i]][[j]][11, ] > 0)/M,
                          mean(R2[[i]][[j]][10, ]), sd(R2[[i]][[j]][10, ]), mean(R2[[i]][[j]][12, ]), sum(R2[[i]][[j]][10, ] - q * R2[[i]][[j]][12, ] < theta2 & R2[[i]][[j]][10, ] + q * R2[[i]][[j]][12, ] > theta2)/M, sum(abs(R2[[i]][[j]][10, ]) - q * R2[[i]][[j]][12, ] > 0)/M,
                          
                          rowMeans(R2[[i]][[j]][13:16, ]))
  }
}

theta1 = 0.2
theta2 = 0.2
T3 = matrix(nrow = 45, ncol = 36)
for (i in 1:5){
  for (j in 1:9){
    T3[(9*(i-1)+j), ] = c(rho[i], err[j],
                          mean(R3[[i]][[j]][1, ]), sd(R3[[i]][[j]][1, ]), mean(R3[[i]][[j]][3, ]), sum(R3[[i]][[j]][1, ] - q * R3[[i]][[j]][3, ] < theta1 & R3[[i]][[j]][1, ] + q * R3[[i]][[j]][3, ] > theta1)/M, sum(abs(R3[[i]][[j]][1, ]) - q * R3[[i]][[j]][3, ] > 0)/M,
                          mean(R3[[i]][[j]][2, ]), sd(R3[[i]][[j]][2, ]), mean(R3[[i]][[j]][4, ]), sum(R3[[i]][[j]][2, ] - q * R3[[i]][[j]][4, ] < theta2 & R3[[i]][[j]][2, ] + q * R3[[i]][[j]][4, ] > theta2)/M, sum(abs(R3[[i]][[j]][2, ]) - q * R3[[i]][[j]][4, ] > 0)/M,
                          
                          mean(R3[[i]][[j]][5, ]), sd(R3[[i]][[j]][5, ]), mean(R3[[i]][[j]][7, ]), sum(R3[[i]][[j]][5, ] - q * R3[[i]][[j]][7, ] < theta1 & R3[[i]][[j]][5, ] + q * R3[[i]][[j]][7, ] > theta1)/M, sum(abs(R3[[i]][[j]][5, ]) - q * R3[[i]][[j]][7, ] > 0)/M,
                          mean(R3[[i]][[j]][6, ]), sd(R3[[i]][[j]][6, ]), mean(R3[[i]][[j]][8, ]), sum(R3[[i]][[j]][6, ] - q * R3[[i]][[j]][8, ] < theta2 & R3[[i]][[j]][6, ] + q * R3[[i]][[j]][8, ] > theta2)/M, sum(abs(R3[[i]][[j]][6, ]) - q * R3[[i]][[j]][8, ] > 0)/M,
                          
                          mean(R3[[i]][[j]][9, ]), sd(R3[[i]][[j]][9, ]), mean(R3[[i]][[j]][11, ]), sum(R3[[i]][[j]][8, ] - q * R3[[i]][[j]][11, ] < theta1 & R3[[i]][[j]][9, ] + q * R3[[i]][[j]][11, ] > theta1)/M, sum(abs(R3[[i]][[j]][9, ]) - q * R3[[i]][[j]][11, ] > 0)/M,
                          mean(R3[[i]][[j]][10, ]), sd(R3[[i]][[j]][10, ]), mean(R3[[i]][[j]][12, ]), sum(R3[[i]][[j]][10, ] - q * R3[[i]][[j]][12, ] < theta2 & R3[[i]][[j]][10, ] + q * R3[[i]][[j]][12, ] > theta2)/M, sum(abs(R3[[i]][[j]][10, ]) - q * R3[[i]][[j]][12, ] > 0)/M,
                          
                          rowMeans(R3[[i]][[j]][13:16, ]))
  }
}

########################################################################################################################
P1 = matrix(nrow = 36, ncol = 8)
for (i in 1:4){
  for (j in 1:9){
    P1[(9*(i-1)+j), ] = c(rho[i], err[j], median(R1[[i]][[j]][1, ]), median(R1[[i]][[j]][2, ]),
                          median(R1[[i]][[j]][5, ]), median(R1[[i]][[j]][6, ]),
                          median(R1[[i]][[j]][9, ]), median(R1[[i]][[j]][10, ]))
  }
}
P2 = matrix(nrow = 36, ncol = 8)
for (i in 1:4){
  for (j in 1:9){
    P2[(9*(i-1)+j), ] = c(rho[i], err[j], median(R2[[i]][[j]][1, ]), median(R2[[i]][[j]][2, ]),
                          median(R2[[i]][[j]][5, ]), median(R2[[i]][[j]][6, ]),
                          median(R2[[i]][[j]][9, ]), median(R2[[i]][[j]][10, ]))
  }
}
P3 = matrix(nrow = 36, ncol = 8)
for (i in 1:4){
  for (j in 1:9){
    P3[(9*(i-1)+j), ] = c(rho[i], err[j], median(R3[[i]][[j]][1, ]), median(R3[[i]][[j]][2, ]),
                          median(R3[[i]][[j]][5, ]), median(R3[[i]][[j]][6, ]),
                          median(R3[[i]][[j]][9, ]), median(R3[[i]][[j]][10 , ]))
  }
}

P_df = as.data.frame(rbind(P1, P2, P3))
colnames(P_df) = c("rho", "err", "IVW1", "IVW2", "Lik1", "Lik2", "Lik1cor", "Lik2cor")
P_df$Scenario = c(rep(paste("S1"), 36), rep("S2", 36), rep("S3", 36))
P_df$rho = rep(c(rep("rho == 0", 9), rep("rho == 0.2", 9), rep("rho == 0.4", 9), rep("rho == 0.6", 9)), 3)
cpal = brewer.pal(4, "Dark2")

P_df = pivot_longer(P_df, 3:8, names_to = "Method")
P_df$Method = factor(P_df$Method, levels = c("IVW1", "IVW2", "Lik1", "Lik2", "Lik1cor", "Lik2cor"))
lab = c(expression(theta[1]^{(IVW)}), expression(theta[2]^{(IVW)}),
        expression(theta[1]^{(MLE)}), expression(theta[2]^{(MLE)}),
        expression(theta[1]^{(MLEcor)}), expression(theta[2]^{(MLEcor)}))
ggplot(as.data.frame(P_df), aes(x = err, y = value)) +
  geom_hline(yintercept = 0.2, color = "grey", size = 0.75) +
  geom_hline(yintercept = 0, color = "grey", size = 0.75) +
  geom_line(aes(color = Method, linetype = Method, size = Method)) +
  scale_color_manual(values = rep(cpal[1:2], 3), labels = lab) +
  scale_linetype_manual(values = c("solid", "solid", "dotted", "dotted", "dashed", "dashed"), labels = lab) +
  scale_size_manual(values = rep(0.75, 6), labels = lab) +
  ylab("Median MR Estimate") + xlab(expression(sigma[zeta[1]]^{2})) +
  facet_grid(rho ~ Scenario, labeller = label_parsed) + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text = element_text(size = 7),
        axis.title = element_text(size = 8), legend.text = element_text(size = 8), strip.text = element_text(size = 8))
