### Variance Componenets Estimation Using the Eigen Decomposition

The code provided below estimates genomic heritability in the following random effects model

   yi= ui +ei 
   
where: ei~iid N(0,vE)  and u=[u1,...,un]' ~N(0,G*vG)

G is an nxn matrix describing covariances between the levels of the random effects. The model assumes one data point per level of the random
effect (this can be easily circunvents).


References: [de los Campos et al., 2010]()

```R
# A function that evaluates the log-likelihood
neg2LogLik<-function(logVar,V,y,d,n=length(y)){
  y<-y-mean(y)
  Vy<-crossprod(V,y)
  Vy2<-as.vector(Vy)^2
  varE<-exp(logVar[1])
  varU<-exp(logVar[2])
  lambda<-varU/varE
  dStar<-(d*lambda+1)
  sumLogD<-sum(log(dStar))
  neg2LogLik_1<- ( n*log(varE) + sumLogD )
  neg2LogLik_2<- (sum(Vy2/dStar))/varE
  out<- neg2LogLik_1+neg2LogLik_2
  return(out)
}

```

### Example 1: Profiling the likelihood

```R
  # Simple simulation
   library(BGLR)
   data(mice)
   
   X=scale(mice.X)
   n=nrow(X) ; p=ncol(X)
   
   h2=0.5
   b=rnorm(sd=sqrt(h2/p),n=p)
   
   signal=X%*%b
   error=rnorm(sd=sqrt(1-h2),n=n)
   y=signal+error
  
   
   G=tcrossprod(X)/p
   EVD=eigen(G)
  
   h2Grid=seq(from=.1,to=.8,by=.01)
   loglik=rep(NA,length(h2Grid))
  
   for(i in 1:length(h2Grid)){
    varG=vP*h2Grid[i]
    varE=vP*(1-h2Grid[i])
    print(c(varE,varG))
    loglik[i]<-neg2LogLik(y=y,V=EVD$vectors,d=EVD$values,logVar=log(c(varE,varG)))
    print(i)
   }
   plot(-loglik~h2Grid,type='l'); abline(v=var(signal)/var(y),col='green')
  
```
### Estimation using optim

```R
    fm<-optim(fn=neg2LogLik,V=EVD$vectors,d=EVD$values,y=y,par=log(c(.2,.8)),
                                n=length(y),hessian=FALSE) 
    varEHat=exp(fm$par[1])
    varGHat=exp(fm$par[2])
    abline(v=(varGHat/(varGHat+varEHat)),col=4)
```


###  Example of an averaging 
