rm(list=ls())


#Libraries(install packages if required)

library(CircStats)
library(foreach)
library(doRNG)
library(progressr)
library(doFuture)
library(future)

registerDoFuture()
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("progress")

#Probability Function

f=function(me,ro){
  delta=pi/4
  betdum=ro*complex(real=cos(me),imaginary=sin(me))
  deldum1=complex(real=cos(-delta),imaginary=sin(-delta))
  deldum2=complex(real=cos(delta),imaginary=sin(delta))
  zerdum=1
  duma=(betdum-deldum1)/(1-Conj(betdum)*deldum1)
  dumarg=Arg(duma)%%(2*pi)
  dumb=(betdum-deldum2)/(1-Conj(betdum)*deldum2)
  dumarg2=Arg(dumb)%%(2*pi)
  dumzer=(betdum-zerdum)/(1-Conj(betdum)*zerdum)
  dumargz=Arg(dumzer)%%(2*pi)
  mi=min(dumarg,dumarg2)
  ma=max(dumarg,dumarg2)
  if(dumargz > mi & dumargz < ma){l1=1-(ma-mi)/(2*pi)}
  else{l1=1-(mi+2*pi-ma)/(2*pi)}
  return(l1)
}

#CARA PROCEDURE

allocation<-function(iniala=10,rtrue,rho=0.92,
                     rho1=0.8,N=200){
  
  trt=c(rep(1,iniala),rep(2,iniala),rep(3,iniala))
  
  if(cov.type=="wrapped"){
    tx1=rwrpcauchy(iniala,0,rho1)
    tx2=rwrpcauchy(iniala,0,rho1)
    tx3=rwrpcauchy(iniala,0,rho1)
  } else {
    tx1=runif(iniala,0,2*pi)
    tx2=runif(iniala,0,2*pi)
    tx3=runif(iniala,0,2*pi)
  }
  
  tx=c(tx1,tx2,tx3)
  ty=mat.or.vec(length(tx),1)
  x<-complex(real=cos(tx),imaginary=sin(tx))
  r<-rtrue[trt]
  tmu<-rwrpcauchy(length(ty),0,rho)
  ty=(Arg((x+r)/(1+r*x))+tmu)%%(2*pi)
  
  ss=length(ty)
  n=length(ty)
  N1=N2=N3=iniala
  
  while(n < N){
    
    fn=function(d){
      rhoe=d[4]
      ret=0
      for(i in 1:ss){
        re=d[trt[i]]
        cmue<-(cos(tx[i])+2*re+re^2*cos(tx[i]))/(1+re^2+2*re*cos(tx[i]))
        smue=-sin(tx[i])*(re^2-1)/(1+re^2+2*re*cos(tx[i]))
        ret=ret+log(1-rhoe^2)-log(1+rhoe^2-
                                    2*rhoe*cos(ty[i])*cmue-2*rhoe*sin(ty[i])*smue)
      }
      return(-ret)
    }
    
    conv<-52; cou=0
    while(conv>0){
      cou=cou+1
      if(cou==1) di <- runif(4,0.01,0.95) else di <- op$par
      
      op=optim(di,fn,
               lower=c(-0.99,-0.99,-0.99,0.01),
               upper=c(0.99,0.99,0.99,0.99),
               method="L-BFGS-B",hessian=(n>50))
      conv<-op$convergence
      if(cou==15){conv=0}
    }
    
    if(cov.type=="wrapped"){
      txn=rwrpcauchy(1,0,rho1)
    } else {
      txn=runif(1,0,2*pi)
    }
    
    xn=complex(real=cos(txn),imaginary=sin(txn))
    me1=Arg((xn+op$par[1])/(1+op$par[1]*xn))
    me2=Arg((xn+op$par[2])/(1+op$par[2]*xn))
    me3=Arg((xn+op$par[3])/(1+op$par[3]*xn))
    
    pme1=f(me1,op$par[4])
    pme2=f(me2,op$par[4])
    pme3=f(me3,op$par[4])
    
    if(n > 50){
      invH = solve(op$hessian)
    } else {
      invH = diag(4)
    }
    
    p1=sqrt(1/pme1)*invH[1,1]*N1
    p2=sqrt(1/pme2)*invH[2,2]*N2
    p3=sqrt(1/pme3)*invH[3,3]*N3
    
    su=p1+p2+p3
    p1=p1/su; p2=p2/su; p3=p3/su
    
    u=runif(1)
    if(u<p1){trtn=1;N1=N1+1}
    else if(u<(p1+p2)){trtn=2;N2=N2+1}
    else{trtn=3;N3=N3+1}
    
    trt=c(trt,trtn)
    tmun=rwrpcauchy(1,0,rho)
    tyn=(Arg((xn+rtrue[trtn])/(1+rtrue[trtn]*xn))+tmun)%%(2*pi)
    
    tx=c(tx,txn)
    ty=c(ty,tyn)
    n=n+1
    ss=n
  }
  
  return(cbind(tx,trt,ty))
}

#CR PROCEDURE

allocation.CR<-function(iniala=10,rtrue,rho=0.92,
                        rho1=0.8,N=200){
  
  trt=c(rep(1,iniala),rep(2,iniala),rep(3,iniala))
  
  if(cov.type=="wrapped"){
    tx1=rwrpcauchy(iniala,0,rho1)
    tx2=rwrpcauchy(iniala,0,rho1)
    tx3=rwrpcauchy(iniala,0,rho1)
  } else {
    tx1=runif(iniala,0,2*pi)
    tx2=runif(iniala,0,2*pi)
    tx3=runif(iniala,0,2*pi)
  }
  
  tx=c(tx1,tx2,tx3)
  ty=mat.or.vec(length(tx),1)
  x<-complex(real=cos(tx),imaginary=sin(tx))
  r<-rtrue[trt]
  tmu<-rwrpcauchy(length(ty),0,rho)
  ty=(Arg((x+r)/(1+r*x))+tmu)%%(2*pi)
  
  N1=N2=N3=iniala
  n=length(ty)
  
  while(n < N){
    
    if(cov.type=="wrapped"){
      txn=rwrpcauchy(1,0,rho1)
    } else {
      txn=runif(1,0,2*pi)
    }
    
    u=runif(1)
    if(u<1/3){trtn=1;N1=N1+1}
    else if(u<2/3){trtn=2;N2=N2+1}
    else{trtn=3;N3=N3+1}
    
    xn=complex(real=cos(txn),imaginary=sin(txn))
    tmun=rwrpcauchy(1,0,rho)
    tyn=(Arg((xn+rtrue[trtn])/(1+rtrue[trtn]*xn))+tmun)%%(2*pi)
    
    trt=c(trt,trtn)
    tx=c(tx,txn)
    ty=c(ty,tyn)
    n=n+1
  }
  
  return(cbind(tx,trt,ty))
}

#LRT

LRT <- function(data){
  
  tx=data[,1]; trt=data[,2]; ty=data[,3]
  ss=length(ty)
  
  fn=function(d){
    rhoe=d[4]; ret=0
    for(i in 1:ss){
      re=d[trt[i]]
      cmue<-(cos(tx[i])+2*re+re^2*cos(tx[i]))/(1+re^2+2*re*cos(tx[i]))
      smue=-sin(tx[i])*(re^2-1)/(1+re^2+2*re*cos(tx[i]))
      ret=ret+log(1-rhoe^2)-log(1+rhoe^2-
                                  2*rhoe*cos(ty[i])*cmue-2*rhoe*sin(ty[i])*smue)
    }
    return(-ret)
  }
  
  di=runif(4,0.01,0.95)
  op=optim(di,fn,
           lower=c(-0.99,-0.99,-0.99,0.01),
           upper=c(0.99,0.99,0.99,0.99),
           method="L-BFGS-B")
  
  L1=-op$value
  
  fn0=function(d){
    rhoe=d[2]; ret=0
    for(i in 1:ss){
      re=d[1]
      cmue<-(cos(tx[i])+2*re+re^2*cos(tx[i]))/(1+re^2+2*re*cos(tx[i]))
      smue=-sin(tx[i])*(re^2-1)/(1+re^2+2*re*cos(tx[i]))
      ret=ret+log(1-rhoe^2)-log(1+rhoe^2-
                                  2*rhoe*cos(ty[i])*cmue-2*rhoe*sin(ty[i])*smue)
    }
    return(-ret)
  }
  
  di=runif(2,0.01,0.95)
  op0=optim(di,fn0,
            lower=c(-0.99,0.01),
            upper=c(0.99,0.99),
            method="L-BFGS-B")
  
  L0=-op0$value
  LR=2*(L1-L0)
  pval=1-pchisq(LR,df=2)
  
  return(pval)
}

LRT.CR <- function(data){ return(LRT(data)) }

#Simulation Functions

AP <- function(it=1000,rtrue,rho=0.8,rho1=0.8,N=200){
  
  with_progress({
    
    p <- progressor(steps = it)
    
    prop <- foreach(i=1:it,.combine=rbind,
                    .packages="CircStats",
                    .options.RNG=100) %dorng% {
                      
                      d <- allocation(rtrue=rtrue,rho=rho,rho1=rho1,N=N)
                      tab <- table(factor(d[,2],levels=1:3))
                      
                      p()   # <-- update progress
                      
                      tab/N
                    }
    
    colMeans(prop)
  })
}

AP.CR <- function(it=1000,rtrue,rho=0.8,rho1=0.8,N=200){
  
  with_progress({
    
    p <- progressor(steps = it)
    
    prop <- foreach(i=1:it,.combine=rbind,
                    .packages="CircStats",
                    .options.RNG=200) %dorng% {
                      
                      d <- allocation.CR(rtrue=rtrue,rho=rho,rho1=rho1,N=N)
                      tab <- table(factor(d[,2],levels=1:3))
                      
                      p()
                      
                      tab/N
                    }
    
    colMeans(prop)
  })
}

pow <- function(it=1000,rtrue,rho=0.8,rho1=0.8,N=200){
  
  with_progress({
    
    p <- progressor(steps = it)
    
    re <- foreach(i=1:it,.combine=c,
                  .packages="CircStats",
                  .options.RNG=300) %dorng% {
                    
                    d <- allocation(rtrue=rtrue,rho=rho,rho1=rho1,N=N)
                    pv <- LRT(d)
                    
                    p()   # update
                    
                    pv < 0.05
                  }
    
    mean(re)
  })
}

pow.CR <- function(it=1000,rtrue,rho=0.8,rho1=0.8,N=200){
  
  with_progress({
    
    p <- progressor(steps = it)
    
    re <- foreach(i=1:it,.combine=c,
                  .packages="CircStats",
                  .options.RNG=400) %dorng% {
                    
                    d <- allocation.CR(rtrue=rtrue,rho=rho,rho1=rho1,N=N)
                    pv <- LRT.CR(d)
                    
                    p()
                    
                    pv < 0.05
                  }
    
    mean(re)
  })
}

error <- function(it=1000,rtrue0=0,rho=0.8,rho1=0.8,N=200){
  r0=c(rtrue0,rtrue0,rtrue0)
  pow(it,r0,rho,rho1,N)
}

error.CR <- function(it=1000,rtrue0=0,rho=0.8,rho1=0.8,N=200){
  r0=c(rtrue0,rtrue0,rtrue0)
  pow.CR(it,r0,rho,rho1,N)
}

set.seed(12345)
cov.type <- "wrapped"

#beta <- c(0,0,0)
#beta<- c(0,-0.2,0.7)
beta<- c(0,-0.1,-0.2)
#beta<- c(0,-0.2,-0.2)

Ns <- c(100,200)
rhos <- c(0.75,0.80)
it <- 1000

results_wrapped <- list()

for(Nval in Ns){
  for(rhoval in rhos){
    
    cat("\n============================\n")
    cat("Wrapped Case | N =",Nval,"| rho =",rhoval,"\n")
    cat("============================\n")
    
    alloc_tab <- rbind(
      CARA = AP(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval),
      CR   = AP.CR(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval)
    )
    
    test_tab <- matrix(
      c(
        pow(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval),
        pow.CR(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval),
        error(it=it,rtrue0=0,rho=rhoval,rho1=0.8,N=Nval),
        error.CR(it=it,rtrue0=0,rho=rhoval,rho1=0.8,N=Nval)
      ),
      nrow=2,byrow=TRUE
    )
    
    rownames(test_tab)<-c("Power","Size")
    colnames(test_tab)<-c("CARA","CR")
    
    results_wrapped[[paste("N",Nval,"rho",rhoval)]] <- list(
      Allocation=alloc_tab,
      Test=test_tab
    )
    
  }
}

cov.type <- "uniform"

results_uniform <- list()

for(Nval in Ns){
  for(rhoval in rhos){
    
    cat("\n============================\n")
    cat("Uniform Case | N =",Nval,"| rho =",rhoval,"\n")
    cat("============================\n")
    
    alloc_tab <- rbind(
      CARA = AP(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval),
      CR   = AP.CR(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval)
    )
    
    test_tab <- matrix(
      c(
        pow(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval),
        pow.CR(it=it,rtrue=beta,rho=rhoval,rho1=0.8,N=Nval),
        error(it=it,rtrue0=0,rho=rhoval,rho1=0.8,N=Nval),
        error.CR(it=it,rtrue0=0,rho=rhoval,rho1=0.8,N=Nval)
      ),
      nrow=2,byrow=TRUE
    )
    
    rownames(test_tab)<-c("Power","Size")
    colnames(test_tab)<-c("CARA","CR")
    
    results_uniform[[paste("N",Nval,"rho",rhoval)]] <- list(
      Allocation=alloc_tab,
      Test=test_tab
    )
    
  }
}

print(results_wrapped)
print(results_uniform)
