library(SuppDists)

rWald <- function(n,B,v,A,s=1)
  # random function for single acumulator
{
  
  rwaldt <- function(n,k,l,s=1,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss
    
    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }
    
    flag <- l>tiny
    x <- rep(NA,times=n)
    
    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- (k/s)^2
    
    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]
    
    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)
    
    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }
  
  # Act as if negative v never terminates, cluge to do single accumulator
  # case by passing negative v
  if (length(v)!=n) v <- rep(v,length.out=n)
  if (length(B)!=n) B <- rep(B,length.out=n)
  if (length(A)!=n) A <- rep(A,length.out=n)
  if (length(s)!=n) s <- rep(s,length.out=n)
  
  # Kluge to return -Inf for negative rates, so can implment one accumulator case
  out <- numeric(n)
  ok <- v>0
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok],s=s[ok])
  out[!ok] <- Inf
  out
}

# dWaldSuppDists <- function(t, k, l, alpha=1) {
#   ## transform k (criterion) / l (rate) parametrisation to nu/lambda parametrisation
#   nu = k / l
#   lambda = (k^2)/(alpha^2)
#   dinvGauss(t, nu=nu, lambda=lambda)
# }
# 

dWald <- function(t,v,B,A,s=1,useSuppDists=TRUE)
  # density for single accumulator
{
  
  digt <- function(t,k=1,l=1,a=.1,s=1,tiny=1e-10) {
    # pdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns digt.0 if a<1e-10
    
    digt.0 <- function(t,k=1,l=1) {
      # pdf of inverse gaussian at t with no k variability
      # much faster than statmod's dinvgauss funciton
      
      lambda <- (k/s)^2
      l0 <- l==0
      e <- numeric(length(t))
      if ( any(!l0) ) {
        mu <- k[!l0]/l[!l0]
        e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
      }
      if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
      x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
      x[t<=0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    if(length(s)!=length(t)) s <- rep(s,length.out=length(t))
    
    tpos <- t <=0
    
    atiny <- a <= tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if ( any(atiny) )
      if(useSuppDists) {
        nu <- k[atiny]/l[atiny]
        lambda <- (k[atiny]/s[atiny])^2
        nu.inf <- is.infinite(nu) | is.na(nu)
        x[atiny][nu.inf] <- 0
        x[atiny][!nu.inf] <- dinvGauss(t[atiny][!nu.inf], 
                                       nu=nu[!nu.inf], 
                                       lambda=lambda[!nu.inf])
      } else {
        x[atiny] <- digt.0(t=t[atiny],k=k[atiny],l=l[atiny],s=s[atiny])
      }
    
    # Threshold variability. CANNOT DO TRIAL-BY-TRIAL s!!
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])
        
        term.2a <- log(.5)+log(l[notltiny])
        term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
        term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.2d <- term.2b+term.2c
        term.2 <- exp(term.2a)*term.2d
        
        term.3 <- term.1+term.2
        term.4 <- log(term.3)-log(2)-log(a[notltiny])
        x[notltiny] <- exp(term.4)
      }
      
      if ( any(ltiny) ) {  # rate zero
        log.t <- log(t[ltiny])
        term.1 <- -.5*(log(2)+log(pi)+log.t)
        term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
        term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
        term.4 <- (exp(-term.2)-exp(-term.3))
        term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
        x[ltiny] <- exp(term.5)
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v>0
  B <- unlist(B)
  A <- unlist(A)
  out[ok] <- digt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2,s=s[ok])
  out[!ok] <- 0
  out
}

pWald <- function(t,v,B,A,s=1, useSuppDists = F)
  # cumulative density for single accumulator
{
  pigt <- function(t,k=1,l=1,a=.1,s=1,tiny=1e-10) {
    # cdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns pigt.0 if a<=0
    
    pigt.0 <- function(t,k=1,l=1,s=1) {
      # cdf of inverse gaussian at t with no k variability
      # much faster than statmod's pinvgauss funciton
      
      mu <- k/l
      lambda <- (k/s)^2
      
      e <- exp(log(2*lambda) - log(mu))
      add <- sqrt(lambda/t) * (1 + t/mu)
      sub <- sqrt(lambda/t) * (1 - t/mu)
      
      p.1 <- 1 - pnorm(add)
      p.2 <- 1 - pnorm(sub)
      x <- exp(e + log(p.1)) + p.2
      
      x[t<0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    if(length(s)!=length(t)) s <- rep(s,length.out=length(t))
    
    tpos <- t<=0
    
    atiny <- a<=tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if ( any(atiny) ){
      if(useSuppDists) {
        nu <- k[atiny]/l[atiny]
        lambda <- (k[atiny]/s[atiny])^2
        nu.inf <- is.infinite(nu) | is.na(nu)
        x[atiny][nu.inf] <- 0
        x[atiny][!nu.inf] <- pinvGauss(t[atiny][!nu.inf], 
                                       nu=nu[!nu.inf], 
                                       lambda=lambda[!nu.inf])
      } else {
        x[atiny] <- pigt.0(t=t[atiny],k=k[atiny],l=l[atiny],s=s[atiny])
      }    
    }
    # Threshold variability
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        log.t <- log(t[notltiny])
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- .5*log.t-.5*log(2*pi)
        term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1 <- exp(term.1a)*(term.1b-term.1c)
        
        term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) +
                         log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) +
                         log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
        
        term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
        term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
        term.4 <- term.4c*term.4a + term.4d*term.4b
        
        x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
      }
      
      if ( any(ltiny) ) {  # rate zero
        sqr.t <- sqrt(t[ltiny])
        log.t <- log(t[ltiny])
        term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
        term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
        term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
        
        term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6 <- 1 + exp(term.6b) - exp(term.6a)
        
        x[ltiny] <- term.5 + term.6
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v>0
  B <- unlist(B)
  A <- unlist(A)
  out[ok] <- pigt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2,s=s[ok])
  out[!ok] <- 0
  out
  
}

rWaldRace <- function(n,v,B,A,t0,s)
  # random function for Wald race.
{
  B[B<0] <- 0 # Protection for negatives
  A[A<0] <- 0
  v[v < 1e-5] <- 1e-5
  bs <- B + runif(length(B), 0, A)
  n_v  <- ncol(v)
  ttf <- matrix(NA, nrow=n, ncol=ncol(v))
  for(i in 1:n_v) {
    ttf[,i] <- rinvGauss(n, nu=bs[,i]/v[,i], lambda=(bs[,i]/s[,i])^2)
  }
  ttf <- ttf + t0
  resp <- apply(ttf, 1, which.min)
  rt <- apply(ttf, 1, min)
  out <- data.frame(RT = rt, R = resp)
  return(out)
}

n1Wald <- function(dt,B,A,t0,s, v, n_v)
  # Generates defective PDF for responses on node=1, dt (decison time) is a vector of times
{
  dt <- dt-t0[,1]  ## assume single t0
  # This makes sure that the Density function is multiplied by the cumulative density function of the other accumululators to get a
  # Defective density function.
  pdf <- dWald(dt, A=A[,1], v=v[,1], B=B[,1], s=s[,1])
  if (n_v > 1) for (i in 2:n_v)
    pdf <- pdf * (1-pWald (dt, A=A[,i], v=v[,i], B=B[,i], s=s[,i]))
  return(pdf)
}

dWaldRace <- function(rt, response, A, B, t0, s, v, st0 = 0, silent = FALSE)
{
  response <- as.numeric(response)
  n_v <- ncol(v)
  out <- vector("numeric", length(rt))
  for (i in unique(response)) {
    sel <- response == i #make sure that only matrix entries belonging to this specific response are selected.
    #Reorders the responses in such a fashion that the current response is the first entry and the other response is the second
    #This reordering is done for all responses.
    out[sel] <- n1Wald(dt = rt[sel],
                       A = A[sel,c(i, seq_len(n_v)[-i])],
                       B = B[sel,c(i, seq_len(n_v)[-i])],
                       t0 = t0[sel,c(i, seq_len(n_v)[-i])],
                       s = s[sel,c(i, seq_len(n_v)[-i])],
                       v = v[sel,c(i, seq_len(n_v)[-i])],
                       n_v= n_v)
  }
  return(out)
}