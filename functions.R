# Functions #

calc <- function(pars) {
  
  m00 <- pars$m00
  m01 <- pars$m01
  m10 <- pars$m10
  m11 <- pars$m11
  
  mu10 <- pars$mu10
  mu00 <- pars$mu00
  mu11 <- pars$mu11
  mu01 <- pars$mu01
  
  
  
  m00t <- function(t) {
    dt <- t-t0
    m00*exp(-1*mu00*dt)/(m00*exp(-1*mu00*dt) + m01*exp(-1*mu01*dt) + m10*exp(-1*mu10*dt) + m11*exp(-1*mu11*dt))
  }
  
  m01t <- function(t) {
    dt <- t-t0
    m01*exp(-1*mu01*dt)/(m00*exp(-1*mu00*dt) + m01*exp(-1*mu01*dt) + m10*exp(-1*mu10*dt) + m11*exp(-1*mu11*dt))
  }
  
  m11t <- function(t) {
    dt <- t-t0
    m11*exp(-1*mu11*dt)/(m00*exp(-1*mu00*dt) + m01*exp(-1*mu01*dt) + m10*exp(-1*mu10*dt) + m11*exp(-1*mu11*dt))
  }
  
  m10t <- function(t) {
    dt <- t-t0
    m10*exp(-1*mu10*dt)/(m00*exp(-1*mu00*dt) + m01*exp(-1*mu01*dt) + m10*exp(-1*mu10*dt) + m11*exp(-1*mu11*dt))
  }
  
  
  t1 <- pars$t1
  t2 <- pars$t2
  t0 <<- t1
  res <- matrix(ncol=11,nrow=0)
  
  k1carr <- 1/(m10 + m11) 
  k1non <- 1/(m01 + m00)
  for(i in t1:t2) {
    dt <- i-t1
    
    k <- m00*exp(-1*mu00*dt) + m01*exp(-1*mu01*dt) + m10*exp(-1*mu10*dt) + m11*exp(-1*mu11*dt)
    m1t <- (m10*exp(-1*mu10*dt) + m11*exp(-1*mu11*dt))/k
    m2t <- (m01*exp(-1*mu01*dt) + m11*exp(-1*mu11*dt))/k
    
    S1carr <- (m10*exp(-1*mu10*dt) + m11*exp(-1*mu11*dt))*k1carr
    S1non <- (m01*exp(-1*mu01*dt) + m00*exp(-1*mu00*dt))*k1non
    
    ld <- round(m11t(i) - (m10t(i) + m11t(i))*(m01t(i) + m11t(i)),8)
    
    p1 <- m10t(i) + m11t(i)
    p2 <- m01t(i) + m11t(i)
    
    q1 <- m01t(i) + m00t(i)
    q2 <- m10t(i) + m00t(i)
    
    r2 <- round(ld^2/(p1*q1*p2*q2),8)
    
    if(r2 < 0) {
      r2 <- -1*r2
    } 
    
    res <- rbind(res, c(i, m00t(i), m01t(i), m11t(i), m10t(i), m1t, m2t, S1carr, S1non, ld, r2))
    
  }
  
  colnames(res) <- c("t", "m00", "m01", "m11", "m10", "m1t", "m2t", "S1carr", "S1non", "ld", "r2")
  
  dd <- list()
  dd$m=res
  dd
}

mu00t <- function(args) {
  mu00 <- args$a*exp(args$b*args$t)
  mu00
}



mu10t <- function(args) {
  if(args$dcase == T) {
    mu10 <- mu00t(args)*(1+args$D1)
  } else {
    mu10 <- mu00t(args)*args$H1
  }
  if(args$epistasis) {
    mu10 <- mu00t(args) + args$R1
  }
  mu10
}

mu01t <- function(args) {
  if(args$dcase == T) {
    mu01 <- mu00t(args)*(1+args$D2)
  } else {
    mu01 <- mu00t(args)*args$H2
  }
  if(args$epistasis) {
    mu01 <- mu00t(args) + args$R2
  }
  mu01
}

mu11t <- function(args) {
  if(args$dcase == T) {
    mu11 <- mu00t(args)*(1+args$D1 + args$D2)
  } else {
    mu11 <- mu00t(args)*args$H1*args$H2
  }
  if(args$epistasis) {
    mu11 <- mu00t(args) + args$R1 + args$R2 + args$c*args$R1*args$R2
  }
  mu11
}


calc_gompertz <- function(pars) {
  print(pars)
  m00 <- pars$m00
  m01 <- pars$m01
  m10 <- pars$m10
  m11 <- pars$m11
  
  a <- pars$a
  b <- pars$b
  
  dcase <- pars$dcase
  D1 <- pars$D1
  D2 <- pars$D2
  H1 <- pars$H1
  H2 <- pars$H2
  
  epistasis <- pars$epistasis
  c <- pars$c
  R1 <- pars$R1
  R2 <- pars$R2
  
  
  m00t <- function(t) {
    pars$t <- t
    dt <- t-t0
    m00*exp(-1*mu00t(pars)*dt)/(m00*exp(-1*mu00t(pars)*dt) + 
                                  m01*exp(-1*mu01t(pars)*dt) + 
                                  m10*exp(-1*mu10t(pars)*dt) + 
                                  m11*exp(-1*mu11t(pars)*dt))
  }
  
  m01t <- function(t) {
    pars$t <- t
    dt <- t-t0
    m01*exp(-1*mu01t(pars)*dt)/(m00*exp(-1*mu00t(pars)*dt) + 
                                  m01*exp(-1*mu01t(pars)*dt) + 
                                  m10*exp(-1*mu10t(pars)*dt) + 
                                  m11*exp(-1*mu11t(pars)*dt))
  }
  
  m11t <- function(t) {
    pars$t <- t
    dt <- t-t0
    m11*exp(-1*mu11t(pars)*dt)/(m00*exp(-1*mu00t(pars)*dt) + 
                                  m01*exp(-1*mu01t(pars)*dt) + 
                                  m10*exp(-1*mu10t(pars)*dt) + 
                                  m11*exp(-1*mu11t(pars)*dt))
  }
  
  m10t <- function(t) {
    pars$t <- t
    dt <- t-t0
    m10*exp(-1*mu10t(pars)*dt)/(m00*exp(-1*mu00t(pars)*dt) + 
                                  m01*exp(-1*mu01t(pars)*dt) + 
                                  m10*exp(-1*mu10t(pars)*dt) + 
                                  m11*exp(-1*mu11t(pars)*dt))
  }
  
  t1 <- pars$t1
  t2 <- pars$t2
  t0 <<- t1
  res <- matrix(ncol=11,nrow=0)
  
  k1carr <- 1/(m10 + m11) 
  k1non <- 1/(m01 + m00)
  
  for(i in t1:t2) {
    dt <- i - t1
    pars$t <- i
    
    k <- m00*exp(-1*mu00t(pars)*dt) + m01*exp(-1*mu01t(pars)*dt) + m10*exp(-1*mu10t(pars)*dt) + m11*exp(-1*mu11t(pars)*dt)
    
    m1t <- (m10*exp(-1*mu10t(pars)*dt) + m11*exp(-1*mu11t(pars)*dt))/k
    
    m2t <- (m01*exp(-1*mu01t(pars)*dt) + m11*exp(-1*mu11t(pars)*dt))/k
    
    S1carr <- (m10*exp(-1*mu10t(pars)*dt) + m11*exp(-1*mu11t(pars)*dt))*k1carr
    
    S1non <- (m01*exp(-1*mu01t(pars)*dt) + m00*exp(-1*mu00t(pars)*dt))*k1non
    
    ld <- round(m11t(i) - (m10t(i) + m11t(i))*(m01t(i) + m11t(i)),8)
    
    p1 <- m10t(i) + m11t(i)
    p2 <- m01t(i) + m11t(i)
    
    q1 <- m01t(i) + m00t(i)
    q2 <- m10t(i) + m00t(i)
    
    denom <- p1*q1*p2*q2
    
    if(!is.na(denom)) {
      if(denom != 0) {
        r2 <- round(ld^2/denom,8)
      } else {
        r2 <- 0
      }
      if(r2 < 0) {
        r2 <- -1*r2
      }
      res <- rbind(res, c(i, m00t(i), m01t(i), m11t(i), m10t(i), m1t, m2t, S1carr, S1non, ld, r2))
    } else {
      r2 <- 0
      ld <- 0
      res <- rbind(res, c(i, 0, 0, m11t(i), m10t(i), m1t, m2t, S1carr, S1non, ld, r2))
    }
  }
  
  colnames(res) <- c("t", "m00", "m01", "m11", "m10", "m1t", "m2t", "S1carr", "S1non", "ld", "r2")
  
  dd <- list()
  dd$m=res
  
  tt <- t1:t2
  pars$t <- tt
  dd[["mu10"]] <- mu10t(pars)
  dd[["mu11"]] <- mu11t(pars)
  dd[["mu01"]] <- mu01t(pars)
  dd[["mu00"]] <- mu00t(pars)
  
  dd
}

