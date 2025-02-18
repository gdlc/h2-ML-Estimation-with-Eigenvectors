

neg2_REML<-function(y,Xt,vars,d,n,nSuccess=NULL,nFailure=NULL){
    ###
    # y=U'eHat where is eHat=OLS residuals
    # Xt=X'U
    # vars=c(varE,varU)
    # d: eigenvalues
    # n: length(y)
    # objective function is -2*L2 (page 1 in http://www.aps.uoguelph.ca/~lrs/Animalz/lesson11/lesson11.pdf)
    # gustavoc@msu.edu 11/30/2016
    ###

    varE=(vars[1])
    varU=(vars[2])
    w<-d*varU+varE
    neg2_reml_1=sum(log(w))
    Xt<-scale(Xt,center=F,scale=1/sqrt(w))
    C=tcrossprod(Xt)
    neg2_reml_2<-log(det(C))
    neg2_reml_3<-sum((y^2)/w)
    h2=varU/(varU+varE)
    if((!is.null(nFailure))&(!is.null(nSuccess))){
      neg2LogPrior=-2* ( (nSuccess-1)*log(h2)+(nFailure-1)*log(1-h2))
    }else{
      neg2LogPrior=0
    }
    out<-neg2_reml_1+neg2_reml_2+neg2_reml_3+neg2LogPrior
    return(out)
}

fitREML<-function(y,EVD=NULL,K=NULL,n=length(y),X=matrix(nrow=n,ncol=1,1),nSuccess=NULL,nFailure=NULL,
                  computeHessian=FALSE,verbose=FALSE){
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


  tmp<-rep(TRUE,length(EVD$values))#
  tmp<-EVD$values>1e-5
  d=EVD$values[tmp]
  U=EVD$vectors[,tmp]

  eHat<-residuals(lsfit(y=y,x=X,intercept=F))
  y<-crossprod(U,eHat)
  Xt<-crossprod(X,U)
  varY<-sum((eHat^2))/(length(y)-ncol(X))

  fm<-bobyqa(par=varY*rep(.5,2),fn=neg2_REML,lower=rep(1e-5,2)*varY,
                                             upper=rep(1.5,2)*varY,
             y=y,Xt=Xt,d=d,n=n,
             nSuccess=nSuccess,nFailure=nFailure)

 estimates=c(fm$par,fm$par[2]/sum(fm$par))
  cornerSol=((estimates[3] <0.001)|((estimates[3] >0.999)))

  if(cornerSol&verbose){ message(' Corner solution, estimates may not be reliable.')}

  names(estimates)<-c('varE','varU','h2')

  out=list(estimates=estimates,logREML=-fm$fval/2)

  out$null_logREML=-neg2_REML(y=y,Xt=Xt,n=n,d=d,nSuccess=nSuccess,nFailure=nFailure,vars=c(varY,0))/2
  out$deviance=out$null_logREML=2*(out$logREML-out$null_logREML)
  out$pVal=pchisq(out$deviance,df=1,lower.tail=FALSE)

 if(computeHessian){
 if(!cornerSol){
   COV1=solve(hessian(fun=neg2_REML,x=fm$par,y=y,Xt=Xt,d=d,n=n,
             nSuccess=nSuccess,nFailure=nFailure))

  X=matrix(nrow=100000,ncol=2,rnorm(200000))%*%chol(COV1)
  X[,1]=X[,1]+estimates[1]
  X[,2]=X[,2]+estimates[2]


  COV=matrix(nrow=3,ncol=3,NA)
  COV[1:2,1:2]=COV1
  h2=X[,2]/rowSums(X)
  COV[3,3]=var(h2)
  COV[3,1]=COV[1,3]=cov(X[,1],h2)
  COV[3,2]=COV[2,3]=cov(X[,2],h2)

  out$SE=sqrt(diag(COV))
  out$vcov=COV
  colnames(out$vcov)=names(estimates)
  rownames(out$vcov)=names(estimates)
  names(out$SE)=names(estimates)
 }else{
   out$SE=rep(NA,3)
   out$vcov=matrix(nrow=3,ncol=3,NA)
   colnames(out$vcov)=names(estimates)
   rownames(out$vcov)=names(estimates)
   names(out$SE)=names(estimates)
  }
 }
  out$msg=fm$msg
  out$nearCornerSol=cornerSol
  return(out)
}

# Example
if(FALSE){
 library(BGLR)
 data(wheat)
 G=tcrossprod(scale(wheat.X))
 G=G/mean(diag(G))
 EVD=eigen(G)
 y=wheat.Y[,2]
 fmREML=fitREML(y=y,EVD=EVD)
 fmML=fitML(y=y,EVD=EVD)
}

