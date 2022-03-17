################################################################################
#Script to re-produce the simulation study relating to mediation analysis
################################################################################

library(MendelianRandomization)
library(parallel)
library(MVMR)
library(tidyverse)
library(RColorBrewer)
library(MASS)

N = 20000
p = 40
q = qnorm(0.975, 0, 1)

M = 1000

cl = makeCluster(6)
clusterEvalQ(cl, library(MendelianRandomization))
clusterEvalQ(cl, library(MVMR))
clusterExport(cl, c('N', 'p', 'M', 'sstat', 'fcalc'))

theta1 = 0
theta2 = 0.2
rho = c(0.6)
err = seq(0, 4, length = 9)
clusterExport(cl, c('theta1', 'theta2', 'rho'))

#X1 null and precise, X2 causal and with error
clusterSetRNGStream(cl, 20220113)
D1_med = lapply(1:9, function(j){
    rr = rho
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
      x2 = rr * x1 + sqrt(1 - rr^2) * (G[, 11:40] %*% bx2[11:40] + e2)
      zeta2 = ee * rnorm(2 * N, 0, 1)
      X1 = x1
      X2 = x2 + zeta2
      Y = theta1 * x1 + theta2 * x2 + U + rnorm(2 * N, 0, 1)
      bxhat1 = bxhat2 = byhat = bzetahat1 = bzetahat2 = vector(length = p)
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

#X2 causal and with error, X1 causal and precise
theta1 = 0.1
theta2 = 0.2
clusterExport(cl, c('theta1', 'theta2'))
clusterSetRNGStream(cl, 20220114)
D2_med = lapply(1:9, function(j){
    rr = rho
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
      x2 = rr * x1 + sqrt(1 - rr^2) * (G[, 11:40] %*% bx2[11:40] + e2)
      zeta2 = ee * rnorm(2 * N, 0, 1)
      X1 = x1
      X2 = x2 + zeta2
      Y = theta1 * x1 + theta2 * x2 + U + rnorm(2 * N, 0, 1)
      bxhat1 = bxhat2 = byhat = bzetahat1 = bzetahat2 = vector(length = p)
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

#X2 causal and with error, X1 null and with error
theta1 = 0
theta2 = 0.2
clusterExport(cl, c('theta1', 'theta2'))
clusterSetRNGStream(cl, 20220115)
D3_med = lapply(1:9, function(j){
    rr = rho
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
      x2 = rr * x1 + sqrt(1 - rr^2) * (G[, 11:40] %*% bx2[11:40] + e2)
      zeta1 = rnorm(2 * N, 0, 1)
      zeta2 = ee * rnorm(2 * N, 0, 1)
      X1 = x1 + zeta1
      X2 = x2 + zeta2
      Y = theta1 * x1 + theta2 * x2 + U + rnorm(2 * N, 0, 1)
      bxhat1 = bxhat2 = byhat = bzetahat1 = bzetahat2 = vector(length = p)
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

#X2 causal and with error, X1 causal and with error
theta1 = 0.1
theta2 = 0.2
clusterExport(cl, c('theta1', 'theta2'))
clusterSetRNGStream(cl, 20220116)
D4_med = lapply(1:9, function(j){
  rr = rho
  ee = sqrt(err[j])
  parLapply(cl, 1:M, function(k){
    g = 0.2
    bx1 = runif(p, 0.08, 0.2)
    bx2 = runif(p, 0.08, 0.2)
    U = rnorm(2 * N, 0, 1)
    maf = runif(p, 0.01, 0.5)
    G = sapply(1:p, function(l){rbinom(2 * N, 2, maf[j])})
    e1 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
    e2 = g * U + sqrt(1 - g^2) * rnorm(2 * N, 0, 1)
    x1 = G %*% bx1 + e1
    x2 = rr * x1 + sqrt(1 - rr^2) * (G[, 11:40] %*% bx2[11:40] + e2)
    zeta1 = rnorm(2 * N, 0, 1)
    zeta2 = ee * rnorm(2 * N, 0, 1)
    X1 = x1 + zeta1
    X2 = x2 + zeta2
    Y = theta1 * x1 + theta2 * x2 + U + rnorm(2 * N, 0, 1)
    bxhat1 = bxhat2 = byhat = bzetahat1 = bzetahat2 = vector(length = p)
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

stopCluster(cl)

set.seed(20220131)
#Estimates on D1_med
R1_med = lapply(1:9, function(j){
    B = sapply(1:M, function(k){
      D = D1_med[[j]][[k]]
      mrob_univ = mr_input(bx = D$bxhat1[1:10], bxse = D$sebx1[1:10], by = D$byhat[1:10], byse = D$seby[1:10])
      mrest_univ = mr_ivw(mrob_univ)
      mrest_me_univ = mrest_me(mrob_univ)
      mrob = mr_mvinput(bx = cbind(D$bxhat1, D$bxhat2), bxse = cbind(D$sebx1, D$sebx2), by = D$byhat, byse = D$seby)
      mrest = mr_mvivw(mrob)
      mrest_me = mrest_me(mrob, no_ini = 1)
      mrest_me_cor = mrest_me_cor(mrob, D$corX, no_ini = 1)
      c(mrest_univ$Estimate, mrest_univ$StdError, mrest_me_univ$thest, sqrt(mrest_me_univ$Var[1, 1]),
        mrest$Estimate, mrest$StdError, mrest_me$thest, sqrt(mrest_me$Var[1, 1]), sqrt(mrest_me$Var[2, 2]),
        mrest_me_cor$thest, sqrt(mrest_me_cor$Var[1, 1]), sqrt(mrest_me_cor$Var[2, 2]), D$fstat)
    })
  })

#Estimates on D2_med
R2_med = lapply(1:9, function(j){
    B = sapply(1:M, function(k){
      D = D2_med[[j]][[k]]
      mrob_univ = mr_input(bx = D$bxhat1[1:10], bxse = D$sebx1[1:10], by = D$byhat[1:10], byse = D$seby[1:10])
      mrest_univ = mr_ivw(mrob_univ)
      mrest_me_univ = mrest_me(mrob_univ)
      mrob = mr_mvinput(bx = cbind(D$bxhat1, D$bxhat2), bxse = cbind(D$sebx1, D$sebx2), by = D$byhat, byse = D$seby)
      mrest = mr_mvivw(mrob)
      mrest_me = mrest_me(mrob, no_ini = 1)
      mrest_me_cor = mrest_me_cor(mrob, D$corX, no_ini = 1)
      c(mrest_univ$Estimate, mrest_univ$StdError, mrest_me_univ$thest, sqrt(mrest_me_univ$Var[1, 1]),
        mrest$Estimate, mrest$StdError, mrest_me$thest, sqrt(mrest_me$Var[1, 1]), sqrt(mrest_me$Var[2, 2]),
        mrest_me_cor$thest, sqrt(mrest_me_cor$Var[1, 1]), sqrt(mrest_me_cor$Var[2, 2]), D$fstat)
    })
  })

#Estimates on D3_med
R3_med = lapply(1:9, function(j){
    B = sapply(1:M, function(k){
      D = D3_med[[j]][[k]]
      mrob_univ = mr_input(bx = D$bxhat1[1:10], bxse = D$sebx1[1:10], by = D$byhat[1:10], byse = D$seby[1:10])
      mrest_univ = mr_ivw(mrob_univ)
      mrest_me_univ = mrest_me(mrob_univ)
      mrob = mr_mvinput(bx = cbind(D$bxhat1, D$bxhat2), bxse = cbind(D$sebx1, D$sebx2), by = D$byhat, byse = D$seby)
      mrest = mr_mvivw(mrob)
      mrest_me = mrest_me(mrob, no_ini = 1)
      mrest_me_cor = mrest_me_cor(mrob, D$corX, no_ini = 1)
      c(mrest_univ$Estimate, mrest_univ$StdError, mrest_me_univ$thest, sqrt(mrest_me_univ$Var[1, 1]),
        mrest$Estimate, mrest$StdError, mrest_me$thest, sqrt(mrest_me$Var[1, 1]), sqrt(mrest_me$Var[2, 2]),
        mrest_me_cor$thest, sqrt(mrest_me_cor$Var[1, 1]), sqrt(mrest_me_cor$Var[2, 2]), D$fstat)
    })
  })

#Estimates on D4_med
R4_med = lapply(1:9, function(j){
  B = sapply(1:M, function(k){
    D = D4_med[[j]][[k]]
    mrob_univ = mr_input(bx = D$bxhat1[1:10], bxse = D$sebx1[1:10], by = D$byhat[1:10], byse = D$seby[1:10])
    mrest_univ = mr_ivw(mrob_univ)
    mrest_me_univ = mrest_me(mrob_univ)
    mrob = mr_mvinput(bx = cbind(D$bxhat1, D$bxhat2), bxse = cbind(D$sebx1, D$sebx2), by = D$byhat, byse = D$seby)
    mrest = mr_mvivw(mrob)
    mrest_me = mrest_me(mrob, no_ini = 1)
    mrest_me_cor = mrest_me_cor(mrob, D$corX, no_ini = 1)
    c(mrest_univ$Estimate, mrest_univ$StdError, mrest_me_univ$thest, sqrt(mrest_me_univ$Var[1, 1]),
      mrest$Estimate, mrest$StdError, mrest_me$thest, sqrt(mrest_me$Var[1, 1]), sqrt(mrest_me$Var[2, 2]),
      mrest_me_cor$thest, sqrt(mrest_me_cor$Var[1, 1]), sqrt(mrest_me_cor$Var[2, 2]), D$fstat)
  })
})

########################################################################################################################
#Tables
########################################################################################################################
T1 = matrix(nrow = 9, ncol = 3)
for (j in 1:9){
  T1[j, ] = apply(sapply(1:M, function(i){
    c(1 - R1_med[[j]][5, i] / R1_med[[j]][1, i],
      1 - R1_med[[j]][9, i] / R1_med[[j]][3, i],
      1 - R1_med[[j]][13, i] / R1_med[[j]][3, i])
    }), 1, median)
}

T2 = matrix(nrow = 9, ncol = 3)
for (j in 1:9){
  T2[j, ] = apply(sapply(1:M, function(i){
    c(1 - R2_med[[j]][5, i] / R2_med[[j]][1, i],
      1 - R2_med[[j]][9, i] / R2_med[[j]][3, i],
      1 - R2_med[[j]][13, i] / R2_med[[j]][3, i])
  }), 1, median)
}

T3 = matrix(nrow = 9, ncol = 3)
for (j in 1:9){
  T3[j, ] = apply(sapply(1:M, function(i){
    c(1 - R3_med[[j]][5, i] / R3_med[[j]][1, i],
      1 - R3_med[[j]][9, i] / R3_med[[j]][3, i],
      1 - R3_med[[j]][13, i] / R3_med[[j]][3, i])
  }), 1, median)
}

T4 = matrix(nrow = 9, ncol = 3)
for (j in 1:9){
  T4[j, ] = apply(sapply(1:M, function(i){
    c(1 - R4_med[[j]][5, i] / R4_med[[j]][1, i],
      1 - R4_med[[j]][9, i] / R4_med[[j]][3, i],
      1 - R4_med[[j]][13, i] / R4_med[[j]][3, i])
  }), 1, median)
}

########################################################################################################################
T_all = as.data.frame(rbind(T1, T2, T3, T4))
T_all$theta1 = c(rep("theta[1] == 0", 9), rep("theta[1] == 0.1", 9), rep("theta[1] == 0", 9), rep("theta[1] == 0.1", 9))
T_all$zeta1 = c(rep("sigma[zeta[1]]^{2} == 0", 18), rep("sigma[zeta[1]]^{2} == 1", 18))
T_all$zeta2 = seq(0, 4, by = 0.5)
T_all$PM0 = c(rep(1, 9), rep(6/11, 9), rep(1, 9), rep(6/11, 9))
names(T_all)[1:3] = c("IVW", "MLE", "MLE (cor)")
T_all = pivot_longer(T_all, 1:3, names_to = "Method")
cpal = brewer.pal(4, "Dark2")

ggplot(data = T_all, aes(x = zeta2, y = value)) + geom_line(aes(color = Method)) +
  facet_grid(theta1 ~ zeta1, labeller = label_parsed, scales = "free_y") + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text = element_text(size = 7),
        axis.title = element_text(size = 8), legend.text = element_text(size = 8), strip.text = element_text(size = 8)) +
  scale_color_manual(values = cpal[1:3]) + geom_hline(aes(yintercept = PM0), lty = 2, size = 0.5) +
  ylab("Median proportion mediated") + xlab(expression(sigma[zeta[2]]^{2}))
