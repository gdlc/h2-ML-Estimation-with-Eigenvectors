### Variance Componenets Estimation Using the Eigen Decomposition

The code provided below estimates genomic heritability in the following random effects model

   yi= ui +ei 
   
where: ei~iid N(0,vE)  and u=[u1,...,un]' ~N(0,G*vG)

G is an nxn matrix describing covariances between the levels of the random effects. The model assumes one data point per level of the random
effect (this can be easily circunvents).


References: [de los Campos et al., 2010]()

```R
# A function that evaluates the log-likelihood
neg2LogLik<-function(logVar,V,d,y){
  ##
  #   logVar: log(residualVariance, geneticVariance)
  #   V: eigenvectors
  #   d: eigenvalues
  #   y: phenotypes
  # Return: the log-likelihood
  ##
  n<-length(y)
  y<-y-mean(y)
  Vy<-crossprod(V,y)
  VySq<-as.vector(Vy)^2
  varE<-exp(logVar[1])
  varG<-exp(logVar[2])
  lambda<-varG/varE
  dStar<-(d*lambda+1)
  sumLogD<-sum(log(dStar))
  logLik_1<- -0.5*( n*log(varE) + sumLogD )
  logLik_2<- (-0.5*sum(VySq/dStar))/varE
  out<- -sum(logLik_1,logLik_2)*2
  return(out)
}
```R

### Example 1: Profiling the likelihood

```R
  library(BGLR)
  data(mice)
  X=scale(mice.X) ## Genotypes
  n=nrow(X) ; p=ncol(X)
  h2=0.5
  b=rnorm(sd=sqrt(h2/p),n=p)
  signal=X%*%b
  error=rnorm(sd=sqrt(1-h2),n=n)
  y=error+signal
  
  
  G=tcrossprod(X)/p
  
  EVD=eigen(G)
  
  vP=var(y)
  
  h2Grid=seq(from=.01,to=.99,by=.01)
  loglik=rep(NA,length(h2Grid))
  
  for(i in 1:length(h2Grid)){
    varG=vP*h2Grid[i]
    varE=vP*(1-h2Grid[i])
    loglik[i]<-neg2LogLik(y=y,V=EVD$vectors,d=EVD$values,logVar=log(c(varE,varG)))
    print(i)
  }
  plot(-loglik~h2Grid,type='l')
  
```
# A wrapper to fit the model

fitVarComp<-function(y,G=NULL,EVD=NULL,minEigVal=1e-4,logVar=log(0.5*rep(var(y),2))){
  if(is.null(EVD)){
    if(is.null(G)){ 
      stop('Either G or EVD must be provided') 
    }else{
      EVD=eigen(G)
    }
  }
  
  tmp=EVD$values>minEigVal
  EVD$vectors=EVD$vectors[,tmp]
  EVD$values=EVD$values[tmp]
  fm=optim(fn=neg2LogLik,par=logVar,y=y,V=EVD$vectors,d=EVD$values)
  
}

```
