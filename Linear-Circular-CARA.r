rm(list = ls())

#Libraries
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

#Probability Function q = P(cos(Y) < cos(delta))
q_fun <- function(mu, rho){
  delta <- pi/4
  
  f <- function(theta){
    (1 - rho^2)/(2*pi*(1 + rho^2 - 2*rho*cos(theta - mu)))
  }
  
  integrate(f, lower = delta, upper = 2*pi - delta)$value
}

#CARA PROCEDURE (Linear–Circular)
allocation <- function(iniala=10, beta.true,
                       rho=0.8, sigma.x=1, N=200){
  
  trt <- c(rep(1,iniala),rep(2,iniala),rep(3,iniala))
  
  # Linear covariate
  x <- rnorm(length(trt), 0, sigma.x)
  
  # Mean direction (mu = 0 baseline)
  mu <- 2 * atan(beta.true[trt] * x)
  
  # Response
  y <- (mu + rwrpcauchy(length(mu),0,rho)) %% (2*pi)
  
  n <- length(y)
  N1 <- N2 <- N3 <- iniala
  
  while(n < N){
    
    #MLE
    fn <- function(d){
      beta <- d[1:3]
      rhoe <- d[4]
      
      ll <- 0
      for(i in 1:n){
        mu_i <- 2 * atan(beta[trt[i]] * x[i])
        ll <- ll + log(1 - rhoe^2) -
          log(1 + rhoe^2 - 2*rhoe*cos(y[i] - mu_i))
      }
      -ll
    }
    
    op <- optim(runif(4,-0.5,0.5), fn,
                lower=c(-2,-2,-2,0.05),
                upper=c(2,2,2,0.99),
                method="L-BFGS-B", hessian=TRUE)
    
    #New patient covariate
    xn <- rnorm(1,0,sigma.x)
    
    muhat <- 2 * atan(op$par[1:3] * xn)
    
    #Failure probabilities
    q <- sapply(muhat, q_fun, rho = op$par[4])
    q <- pmax(q, 1e-6)
    
    #Variance (Fisher info inverse)
    H <- (op$hessian + t(op$hessian))/2
    
    V <- tryCatch(
      solve(H),
      error=function(e) diag(1,4)
    )
    
    Vdiag <- pmax(diag(V)[1:3], 1e-6)
    
    #CARA allocation rule
    p_raw <- sqrt(Vdiag / q)
    
    if(any(!is.finite(p_raw)) || sum(p_raw) <= 0){
      p <- rep(1/3,3)
    } else {
      p <- p_raw / sum(p_raw)
    }
    
    #Assign treatment
    u <- runif(1)
    trtn <- which.max(u <= cumsum(p))
    
    if(trtn==1) N1 <- N1+1
    if(trtn==2) N2 <- N2+1
    if(trtn==3) N3 <- N3+1
    
    #Generate response
    trt <- c(trt, trtn)
    x   <- c(x, xn)
    
    mu.true <- 2 * atan(beta.true[trtn] * xn)
    y <- c(y, (mu.true + rwrpcauchy(1,0,rho)) %% (2*pi))
    
    n <- n + 1
  }
  
  cbind(x, trt, y)
}

#Complete Randomization
allocation.CR <- function(iniala=10, beta.true,
                          rho=0.8, sigma.x=1, N=200){
  
  trt <- c(rep(1,iniala),rep(2,iniala),rep(3,iniala))
  x <- rnorm(length(trt), 0, sigma.x)
  
  mu <- 2 * atan(beta.true[trt] * x)
  y  <- (mu + rwrpcauchy(length(mu),0,rho)) %% (2*pi)
  
  n <- length(y)
  
  while(n < N){
    
    xn <- rnorm(1,0,sigma.x)
    trtn <- sample(1:3,1)
    
    trt <- c(trt, trtn)
    x   <- c(x, xn)
    
    mu <- 2 * atan(beta.true[trtn] * xn)
    y  <- c(y, (mu + rwrpcauchy(1,0,rho)) %% (2*pi))
    
    n <- n + 1
  }
  
  cbind(x, trt, y)
}

#Likelihood Ratio Test
LRT <- function(data){
  
  x <- data[,1]; trt <- data[,2]; y <- data[,3]
  n <- length(y)
  
  fn0 <- function(d){
    ll <- 0
    for(i in 1:n){
      mu <- 2 * atan(d[1]*x[i])
      ll <- ll + log(pmax(1 - d[2]^2, 1e-8)) -
        log(pmax(1 + d[2]^2 - 2*d[2]*cos(y[i]-mu), 1e-8))
    }
    -ll
  }
  
  fn1 <- function(d){
    ll <- 0
    for(i in 1:n){
      mu <- 2 * atan(d[trt[i]]*x[i])
      ll <- ll + log(pmax(1-d[4]^2,1e-8)) -
        log(pmax(1+d[4]^2-2*d[4]*cos(y[i]-mu),1e-8))
    }
    -ll
  }
  
  op0 <- optim(c(0,0.5), fn0, method="L-BFGS-B",
               lower=c(-2,0.05), upper=c(2,0.99))
  
  op1 <- optim(runif(4), fn1, method="L-BFGS-B",
               lower=c(-2,-2,-2,0.05),
               upper=c(2,2,2,0.99))
  
  LR <- 2*(op0$value - op1$value)
  pval <- 1 - pchisq(LR, df=2)
  
  return(pval)
}

#Simulation Functions

AP <- function(it=500, beta.true, rho=0.8,
               sigma.x=1, N=200){
  
  with_progress({
    p <- progressor(steps = it)
    
    prop <- foreach(i=1:it, .combine=rbind,
                    .packages="CircStats",
                    .options.RNG=100) %dorng% {
                      
                      d <- allocation(beta.true=beta.true,
                                      rho=rho,
                                      sigma.x=sigma.x,
                                      N=N)
                      
                      tab <- table(factor(d[,2],levels=1:3))
                      p()
                      tab/N
                    }
    
    colMeans(prop)
  })
}

AP.CR <- function(it=500, beta.true, rho=0.8,
                  sigma.x=1, N=200){
  
  with_progress({
    p <- progressor(steps = it)
    
    prop <- foreach(i=1:it, .combine=rbind,
                    .packages="CircStats",
                    .options.RNG=200) %dorng% {
                      
                      d <- allocation.CR(beta.true=beta.true,
                                         rho=rho,
                                         sigma.x=sigma.x,
                                         N=N)
                      
                      tab <- table(factor(d[,2],levels=1:3))
                      p()
                      tab/N
                    }
    
    colMeans(prop)
  })
}

pow <- function(it=500, beta.true, rho=0.8,
                sigma.x=1, N=200){
  
  with_progress({
    p <- progressor(steps = it)
    
    re <- foreach(i=1:it, .combine=c,
                  .packages="CircStats",
                  .options.RNG=300) %dorng% {
                    
                    d <- allocation(beta.true=beta.true,
                                    rho=rho,
                                    sigma.x=sigma.x,
                                    N=N)
                    
                    pv <- LRT(d)
                    p()
                    pv < 0.05
                  }
    
    mean(re)
  })
}

pow.CR <- function(it=500, beta.true, rho=0.8,
                   sigma.x=1, N=200){
  
  with_progress({
    p <- progressor(steps = it)
    
    re <- foreach(i=1:it, .combine=c,
                  .packages="CircStats",
                  .options.RNG=400) %dorng% {
                    
                    d <- allocation.CR(beta.true=beta.true,
                                       rho=rho,
                                       sigma.x=sigma.x,
                                       N=N)
                    
                    pv <- LRT(d)
                    p()
                    pv < 0.05
                  }
    
    mean(re)
  })
}

size <- function(it=500, rho=0.8,
                 sigma.x=1, N=200){
  
  pow(it=it,
      beta.true=c(0,0,0),
      rho=rho,
      sigma.x=sigma.x,
      N=N)
}

size.CR <- function(it=500, rho=0.8,
                    sigma.x=1, N=200){
  
  pow.CR(it=it,
         beta.true=c(0,0,0),
         rho=rho,
         sigma.x=sigma.x,
         N=N)
}

#Test Cases
beta <- c(0, -0.1, 0.2)
Ns <- c(100, 200)
rhos <- c(0.75, 0.80)
it <- 200
results_LC <- list()
set.seed(123)
for(Nval in Ns){
  for(rhoval in rhos){
    
    cat("\n============================\n")
    cat("Linear-Circular | N =",Nval,"| rho =",rhoval,"\n")
    cat("============================\n")
    
    #Allocation Table
    
    alloc_tab <- rbind(
      CARA = AP(it=it,
                beta.true=beta,
                rho=rhoval,
                sigma.x=1,
                N=Nval),
      
      CR   = AP.CR(it=it,
                   beta.true=beta,
                   rho=rhoval,
                   sigma.x=1,
                   N=Nval)
    )
    
    #Power & Size Table
    
    test_tab <- matrix(
      c(
        pow(it=it,
            beta.true=beta,
            rho=rhoval,
            sigma.x=1,
            N=Nval),
        
        pow.CR(it=it,
               beta.true=beta,
               rho=rhoval,
               sigma.x=1,
               N=Nval),
        
        size(it=it,
             rho=rhoval,
             sigma.x=1,
             N=Nval),
        
        size.CR(it=it,
                rho=rhoval,
                sigma.x=1,
                N=Nval)
      ),
      nrow=2, byrow=TRUE
    )
    
    rownames(test_tab) <- c("Power","Size")
    colnames(test_tab) <- c("CARA","CR")
    
    #Store results
    
    results_LC[[paste("N=",Nval,"rho=",rhoval)]] <- list(
      Allocation = alloc_tab,
      Test = test_tab
    )
    
  }
}

results_LC <- lapply(results_LC, function(x){
  x$Allocation <- round(x$Allocation,3)
  x$Test <- round(x$Test,3)
  x
})

print(results_LC)

