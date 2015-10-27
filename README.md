### Variance Componenets Estimation Using the Eigen Decomposition
*Contact*: Gustavo de los Campos (gdeloscampos@gmail.com)


The code presented below contains functions for evaluating the likelihood of a random effects model using the eigenvalue decomposition of the co-variance matrix associated to a random effect. The model and derivations are described briefly [here](https://github.com/gdlc/h2-ML-Estimation-with-Eigenvectors/blob/master/simple_neg2loglik.pdf). 

Several authors have explited the equivalence between Gaussian processes and reandom regressions on eigenvectors, a few examples of these are: [de los Campos et al. (2010)](http://www.ncbi.nlm.nih.gov/pubmed/20943010) ,[Zhou and Stephens (2012) ](http://www.ncbi.nlm.nih.gov/pubmed/22706312?dopt=Abstract&holding=npg) and  [Janss et al., (2012)](http://www.genetics.org/content/192/2/693.full.pdf). 


### An R function for evaluating the log-likelihood.
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

### Example 1: Evaluation of the likelihood over a grid of values of heritability
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
   
   vP=var(y)
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

### Two-dimension grid search

```R
    library(graphics)
    vP=var(y)
    vE=seq(from=vP/10,to=vP*9/10,by=.01)
    vG=vE
    
    OUT=matrix(nrow=length(vE),ncol=length(vG),0)
    colnames(OUT)=paste0('vG=',round(vG,4))
    rownames(OUT)=paste0('vE=',round(vE,4))
    
    for(i in 1:length(vE)){
      for(j in 1:length(vG)){
        OUT[i,j]<-neg2LogLik(y=y,V=EVD$vectors,d=EVD$values,logVar=log(c(vE[i],vG[j])))
      }
    }
    contour(x=vE,y=vG,z=OUT,nlevels=400,xlab='vE',ylab='vG')
    points(x=exp(fm$par[1]),y=exp(fm$par[2]),col=2,pch=19,cex=1.5)

    abline(v=(varGHat/(varGHat+varEHat)),col=4)
```

### Likelihood profiling

The functions below can be used to proflie a normal likelihood relative to values of heritability. For each value of heritability in a user-defined grid, the function profile_h2 obtains the maximum likelihood estimator of the phenotypic variance and evaluates the log-likelihood for the value of h2 and the MLE of phenotypic variance. 

#### A function to evaluates the likelihood as a function of phenotypic variance and heritability

```R
## This function evaluates the log-likelihood as a function of the phenotypic variance and heritability
neg2LogLik_h2<-function(varP,h2,V,y,d,n=length(y)){
  # evaluates the log-likeihood as a function of h2 and varP (phentoypic variance)
  y<-y-mean(y)
  Vy<-crossprod(V,y)
  Vy2<-as.vector(Vy)^2
  varE<-varP*(1-h2)
  varU<-varP*h2
  lambda<-varU/varE
  dStar<-(d*lambda+1)
  sumLogD<-sum(log(dStar))
  neg2LogLik_1<- ( n*log(varE) + sumLogD )
  neg2LogLik_2<- (sum(Vy2/dStar))/varE
  out<- neg2LogLik_1+neg2LogLik_2
  return(out)
}
```

#### Use

```R
  library(BGLR)
  data(wheat)
  y<-wheat.Y[,1]
  G<-tcrossprod(scale(wheat.X))
  G<-G/mean(diag(G))
  EVD<-eigen(G)
  V=EVD$vectors[,EVD$values>1e-5]
  d<-EVD$values[EVD$values>1e-5]
  
  neg2LogLik_h2(y=y,V=V,d=d,varP=1.1,h2=.4)

```

#### A funciton for profiling the likelihood

```R
profile_h2<-function(y,V,d,n=length(y),h2=seq(from=1/100,to=I(1-1/100),by=1/1000),plot=TRUE,returnResults=T){
 	varP_int<-c(.5,2)*var(y)
	logLik=rep(NA,length(h2))
	for(i in 1:length(h2)){
		fm=optimize(f=neg2LogLik_h2,V=EVD$vectors,d=EVD$values,y=y,h2=h2[i],
                                n=length(y),interval=varP_int) 
        logLik[i]= -2*fm$objective
	}
	
	tmp<-which(logLik==max(logLik))
	cond1<-logLik>max(logLik)-10
	cond2<-logLik<max(logLik)+10
	tmp<-(cond1&cond2)
	x2<-logLik[tmp]
	x1<-h2[tmp]
  if(plot){  
    plot(x2~x1,xlab='h2',ylab='Log-Likelihood',type='l',col=4)
    abline(v=h2[which(logLik==max(logLik))],col=2)
  }
  if(returnResults){ return( data.frame(h2=h2,logLik=logLik)) }
}
```

#### Use

```R
  tmp=profile_h2(y=y,V=V,d=d)
  head(tmp)
```
