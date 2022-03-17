################################################################################
#Functions used in the simulation studies
################################################################################

fcalc = function(X1, X2, G){
  f1 = summary(lm(X1 ~ G))$f[1]
  f2 = summary(lm(X2 ~ G))$f[1]
  x2fit = lm(X2 ~ G)$fitted
  x1resid = lm(X1 ~ x2fit)$resid
  cf1 = summary(lm(x1resid ~ G))$f[1]
  x1fit = lm(X1 ~ G)$fitted
  x2resid = lm(X2 ~ x1fit)$resid
  cf2 = summary(lm(x2resid ~ G))$f[1]
  c(f1, f2, cf1, cf2)
}

mv_norm = function(n, mu, Sigma){
  d = dim(Sigma)[1]
  if (length(mu) != d){
    stop('mu and Sigma must be the same dimension.')
  }
  A = chol(Sigma)
  Z = sapply(1:d, function(j){rnorm(n, 0, 1)})
  M = matrix(rep(mu, n), nrow = n, byrow = TRUE)
  X = drop(M) + Z %*% A
}

sstat = function(Y, X, intercept = TRUE){
  n = length(Y)
  if (intercept == TRUE){
    xx = cbind(rep(1, n), X)
  }
  else {xx = X}
  mod = lm.fit(xx, Y)
  bhat= c(mod$coefficients[2])
  s = t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)
  se = sqrt((c(s) * solve(t(xx) %*% xx))[2,2])
  return(list("bhat" = bhat, "se" = se))
}