neg2_REML<-function(eStar,Xt,vars,d,n){
	###
    # eStar=V'e=V'(y-Xb)
    # Xt=X'V
    # vars=c(varE,varU)
    # d: eigenvalues
    # n: length(y)
    ###     
	varE=vars[1]
	varU=vars[2]
	lambda<-varU/varE
	w<-d*lambda+1
	neg2_reml_1<- n*log(varE)+sum(log(w))
		
	Xt<-scale(Xt,center=F,scale=1/sqrt(w))
	XVInvX=tcrossprod(Xt)
	
	neg2_reml_2<-log(det(XVInvX))
	
	SS<-sum((eStar^2)/w)
	neg2_reml_3<-SS/varE
	
	out<-neg2_reml_1+neg2_reml_2+neg2_reml_3
	return(out)
}

fitREML <- function(y, EVD = NULL, K = NULL, n = length(y), X = matrix(nrow = n, ncol = 1, 1),
                  nSuccess = 1, nFailure = 1, computeHessian = TRUE){
 
  ###
  # y: phenotype
  # EVD=eigen(G)
  # X: incidence matrix for fixed effects
  # gustavoc@msu.edu 11/30/2016
  ###
 
  if(is.null(EVD)){
    if(is.null(K)){
      stop('provide either K or its eigenvalue decomposition')
    }else{
      EVD=eigen(K)
    }
  }
  library(minqa)
  library(numDeriv)
 
  varY <- var(y)
  tmp <- rep(TRUE,length(EVD$values))#EVD$values>1e-5
  d = EVD$values[tmp]
  U = EVD$vectors[,tmp]
  eHat <- residuals(lsfit(y = y, x = X, intercept = F))
  y <- crossprod(U, eHat)
  Xt <- crossprod(X, U)
  fm <- bobyqa(par = varY * rep(.5, 2), fn = neg2_REML, lower = rep(1e-5, 2) * varY, upper = rep(1.5, 1.5) * varY, eStar= y, Xt = Xt, d = d, n = n)
  estimates = c(fm$par, fm$par[2] / sum(fm$par))
  names(estimates) <- c('varE', 'varU', 'h2')
  out = list(estimates = estimates, logREML = -2 * fm$fval)
  if(computeHessian){
    COV1 = solve(hessian(fun = neg2_REML, x = fm$par, eStar = y, Xt = Xt, d = d, n = n))
    X = matrix(nrow = 100000, ncol = 2, rnorm(200000)) %*% chol(COV1)
    X[,1] = X[,1] + estimates[1]
    X[,2] = X[,2] + estimates[2]
    COV = matrix(nrow = 3, ncol = 3, NA)
    COV[1:2, 1:2] = COV1
    h2 = X[,2] / rowSums(X)
    COV[3, 3] = var(h2)
    COV[3, 1] = COV[1, 3] = cov(X[,1], h2)
    COV[3,2] = COV[2,3] = cov(X[,2], h2)
    out$SEs = sqrt(diag(COV))
    out$vcov = COV
  }
  return(out)
}
