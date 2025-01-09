# This code creates functions for simultaneouse Tolerance intervals

###############################################################
###############################################################
# Product set 

# P.target is a target contnet proportion.
# alpha is our 1 - confidence level (gamma). 
# X is our design matrix.
# x0 is a vector of x points we would like to draw observation from.

PS_k_factor= function(P.target,alpha.target,X,x0,adjust=F){
  
P = P.target

if (adjust == F){
  alpha = alpha.target
}  else if (adjust == T){
  alpha = 2*alpha.target
}


n = dim(X)[1]
q = dim(X)[2]

  # k1 is defined in Limam 1988 page802
  k1.func = function(q,n,alpha){
    value = sqrt(q*qf(1-alpha/2,q,n-q))
    return(value)
  }
  
  # k2 is defined in Limam 1988 page802
 
  k2.func = function(q,n,alpha){
    value = sqrt(qchisq(alpha/2,n-q)/(n-q))
    return(value)
  }
  
  # for each x point j, find a k-factor
  
  k.factor = rep(NA,dim(x0)[1])
  
  for (j in 1:length(k.factor)){
    
    # x.point is the j^th row of x0
    x.point = x0[j,]
    
    # delta is defined in Limam 1988 page802
    delta.func = function(x.point){
      value = as.numeric(sqrt(t(x.point)%*%solve(t(X)%*%X)%*%x.point))
      return(value)
    }
    
    # C(k2*k1*delta, k2*k.factor) is defined in Limam 1988 page802
    C.func = function(k.factor,k1,k2,delta){
      value = pnorm(k1*k2*delta+k2*k.factor) - pnorm(k1*k2*delta-k2*k.factor)
      return(value)
    }
    
    max.cand = 50
    k.candidate = seq(0,max.cand,by=0.01)
    C.minus.P.vec = rep(NA,length(k.candidate))
    
    for (i in 1:length(k.candidate)){
      
      k1 = k1.func(q,n,alpha)
      k2 = k2.func(q,n,alpha)
      delta = delta.func(x.point)
      
      k.factor[j]=k.candidate[i]
      C =  C.func(k.factor[j],k1,k2,delta)
      
      C.minus.P.vec[i] = abs(C-P)
      #print(i)
    } # end of for loop i  
    
    k.factor[j]  = k.candidate[which.min(C.minus.P.vec)]
    
  } # end of for loop j
  
  if (max(k.factor) == max.cand){
    warning("k-factor is the max candidate value. Ajudst the max value.")
  }
  
  return(k.factor)
  # output = list(,delta)
  # names(output) = c("kfactors","delta")
  # return(output)
}

###############################################################
###############################################################
# Chevostekova2 method

C2_k_factor= function(P.target,alpha.target,X,x0,h_sq=NULL){
  
  X=X
  x0=x0
  P = P.target
  alpha = alpha.target
  gamma = 1-alpha.target
  
  if (is.null(h_sq)){ 
    
    # If h_sq is not specified, compute it by using procedure
    # introduced in Limam's dissertion.
    
    # a is a grid of values from 0 to 1 by 0.1.
    ai = seq(0,1,0.1)
    
    # r_i is r that makes r_P = 0 for each ai
    r_P = function(a,r,P) pnorm(a+r) - pnorm(a-r) - P
    ri = rep(NA,length(ai))
    
    for (i in 1:length(ai)){
      
      # rp is the value in the candidates r.cand that makes r_P = 0   
      r.cand = seq(0,5,0.001)
      rp = r.cand[which.min(abs(r_P(ai[i],r.cand,P)))]
      ri[i] = rp
      
    } # end  for (i in 1:length(ai))
    
    # h is chosen to mimize the following sum of sqaured deviation (SSD).
    # h2 = h^2
    
    SSD = function(h2,P){
      values = (ri - rep(qnorm(P),length(ai)) - sqrt(rep(h2,length(ai)) + ai^2))^2
      return(sum(values))
    }
    
    SSDV = Vectorize(SSD)
    sum_of_sqd = SSDV(seq(0,1,0.0001),P)
    
    # h2 is the value of h2 that minimizes SSDV 
    h2 = seq(0,1,0.0001)[which.min(sum_of_sqd)]
    
  } else {
    h2 = h_sq
  }
  
  n = dim(X)[1]
  q = dim(X)[2]
  
  find_zp = function(P) qnorm(P)
  find_rm = function(zp,n,d) return((zp + sqrt(zp^2 + 4*n*d^2))/2)
  
  # When n <= 60, use the formula from Chvostekova's paper 
  # to find the gamma quantile.
  
  find_m =function(alpha,X){
    
    n=dim(X)[1]
    q=dim(X)[2]
    
    difference = function(m){
      func = function(w) pchisq(q=m+n*log(w),df=q)*dchisq(x=w,df=n-q)
      G = integrate(func,lower =0, upper = Inf)$value
      return(G-gamma)
    }
    
    V_difference = Vectorize(difference)
    
    answer = uniroot(V_difference,interval = c(-1000,1000))$root
    m_value = round(answer,3)
    # plot(seq(-100,100,by=0.1),abs(V_difference(seq(-100,100,by=0.1))),type="l")
    # abline(v=m_value,col="red")
    return(m_value)
    
  } # end find_m()
  
  M_value = function(alpha,X){
    
    n=dim(X)[1]
    q=dim(X)[2]
    
    values = rchisq(n=100000,df=q) - n*log(rchisq(n=100000,df=n-q))
    #hist(values)
    m_val = quantile(values,probs=1-alpha)
    #abline(v=m_val, col="red")
    return(m_val)
    
  }
  
  # for each x point j, find a k-factor
  
  k.factor = rep(NA,dim(x0)[1])
  
  # Find the value of m, which is the gamma^th quantile of M
  
  if (dim(X)[1] <=60){
    m = find_m(alpha=alpha,X=X)
  } else if (dim(X)[1] > 60){
    mmm = rep(NA,1000)
    for (i in 1:1000){
      mmm[i]=M_value(X=X,alpha=alpha)
      # print(i)
    } # end for (i in 1:1000)
    m = mean(mmm)
  }
  
  
  for (j in 1:length(k.factor)){
    
    # x.point is the j^th row of x0
    x.point = x0[j,]
    
    # delta is defined in Limam 1988 page802
    delta.func = function(x.point){
      value = as.numeric(sqrt(t(x.point)%*%solve(t(X)%*%X)%*%x.point))
      return(value)
    }
    
    d = delta.func(x.point=x.point)
    zp = find_zp(P)
    rm = find_rm(zp,n,d)
    
   # m was used to be computed here, but let's get it before the for loop,
   # since it only dependes on X and alpha.  
    
    a = (rm - zp)^2 - h2 - (d^2)*m 
    b = n*d^2
    
    k.factor[j]  = sqrt(((rm^2)*(n-q))/exp(a/b))
    
  } # end of for loop j
  output = list(round(k.factor,2),h2,m,d)
  names(output) = c("kfactors","h2","m","delta")
  return(output)
}



