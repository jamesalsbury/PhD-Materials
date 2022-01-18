glm1 <- function(y, X, bp, bpp){
  BetaVec <- rep(1, ncol(X))
  bppfunc <- bpp(X%*%BetaVec)
  zk <- y - bp(X%*%BetaVec)
  BetaVec <- BetaVec + solve(t(X)*bppfunc*X)*t(X)*zk
}

poiReg <- function(formula, data){
  mf=model.frame(formula, data = data)
  y=model.response(mf)
  X=model.matrix(formula, mf)
  glm1(y, X, exp, exp)
}


poiReg(Claims ~ log(Holders), data = MASS::Insurance)
glm(Claims ~ log(Holders), data = MASS::Insurance, family = "poisson")


