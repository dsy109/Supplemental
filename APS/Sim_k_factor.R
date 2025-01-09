# This code creates functions for simultaneouse Tolerance intervals
# 1, Product Set approach,
# 2, Wilson's approach,
# 3, Modified Wilson's approach,

# x.values = function(number_of_points){
#   return(cbind(1,seq(from = 1.31, to = 1.4,length.out=number_of_points)))
# }
# 
# design_matrix = function(number_of_points,replication){
#   return(cbind(1,sort(rep(seq(from=1.31,to=1.4,length.out=number_of_points),replication))))
# }
# 
# X = design_matrix(number_of_points=3,replication=1)
# x0 = x.values(number_of_points=3)

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
# Wilson

# P.target = 0.75
# alpha.target = 0.95

Wilson_k_factor = function(P.target, alpha.target,X, x0, h_sq=NULL){
  
  P = P.target
  r0 = qnorm(P)
  
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
  
  alpha = alpha.target
  
  delta.func = function(x.point){
    value = as.numeric(t(x.point)%*%solve(t(X)%*%X)%*%x.point)
    return(value)
  }
  
  k_factor = rep(NA,length(x0[,1]))
  
  for (j in 1:length(x0[,1])){
    
    x.point=x0[j,]
    delta = delta.func(x.point)
    
    n = dim(X)[1]
    q = dim(X)[2]
    v = n-q
    k = sqrt((2*v-1)/(2*v))
    c = qchisq(p=1-alpha,df=q+1) 
    g = (r0^2) - h2 - delta*(c-2*v*k^2)
    
    A = 4*(r0^2) - 4*g
    B = 16*r0*v*k*delta
    C = -8*v*delta*g + 16*(v^2)*(k^2)*(delta^2)
    
    if (B^2 - 4*A*C < 0){　# If the inside of the sqrt is negative
      warning("B^2-$AC is less than 0. It is treated as 0")
      k_factor[j] = round(max(c((-B)/(2*A),(-B)/(2*A))),2)
    } else if (B^2 - 4*A*C >= 0){
      k_factor[j] = round(max(c((-B - sqrt(B^2 - 4*A*C))/(2*A),(-B + sqrt(B^2 - 4*A*C))/(2*A))),2)
    }
    
    # k_factor[j] = round(max(c((-B - sqrt(B^2 - 4*A*C))/(2*A),(-B + sqrt(B^2 - 4*A*C))/(2*A))),2)
    
  } # end of for loop
  
  # Take the absolute value of the calculated k-factors 
  # in case it is negative.
  k_factor = abs(k_factor)
  
  #if(min(k_factor)<0){
  #  warning("There is at least one k factor less than 0")
  #}
  
  output = list(k_factor,h2,c)
  names(output) = c("kfactors","h2","c")
  return(output)
  
}

# X = design_matrix(number_of_points=3,replication=1)
# x0 = x.values(number_of_points=3)
# Wilson_k_factor(0.75,0.95,X,x0)
# 
# X = design_matrix(number_of_points=3,replication=2)
# x0 = x.values(number_of_points=3)
# Wilson_k_factor(0.75,0.95,X,x0)
# 
# X = design_matrix(number_of_points=4,replication=1)
# x0 = x.values(number_of_points=4)
# Wilson_k_factor(0.75,0.95,X,x0)

###############################################################
###############################################################
# Modified Wilson


# cm depends on alpha and X, but it takes time to compute.
# So here, we just put the value of cm

MW_k_factor = function(P.target, alpha.target,X, x0, cm=NULL, h_sq=NULL){
  
  alpha = alpha.target
  n = dim(X)[1]
  q = dim(X)[2]
  v = n-q
  k = sqrt((2*v-1)/(2*v))
  
  if (is.null(cm)){
   
    find_cm2 =function(gamma,X){
    
      difference = function(c){
        cm=c
        func1 = function(w1) pchisq(q=w1*cm/(v*k^2),df=q)*dchisq(x=w1,df=v)
        func2 = function(w2) pchisq(q=cm-2*v*(sqrt(w2/v)-k)^2,df=q)*dchisq(x=w2,df=v)
        G1 = integrate(func1,lower = v*k^2, upper = Inf)$value
        G2 = integrate(func2,lower = v*(k-sqrt(cm/(2*v)))^2, upper = v*k^2)$value
        return(G1+G2-gamma)
      }
      
      differenceV = Vectorize(difference)
      
      # curve(differenceV,xlim=c(0,100),ylim=c(-1,0.2))
      # abline(h=0,col="red")
      
      if (max(differenceV(seq(0,100,by=0.1)))<0){
        warning("There is no c_m that satisfies the equation.")
        cm = seq(0,20,by=0.001)[which.max(differenceV(seq(0,20,by=0.001)))]
        abline(v=cm,col="blue")
        return(cm)
      } else {
        return(round(uniroot(differenceV,interval = c(0,20))$root,3))
      }
      
    } # end find_cm2() function 
    
    # Calculated the value of cm 
    c = find_cm2(gamma=1-alpha,X=X)
    
  } else {
    c = cm
  } # end if (is.null(cm))
  
  P = P.target
  r0 = qnorm(P)
  
  if (is.null(h_sq)){ 
    h2 = round((qnorm((P+1)/2) - qnorm(P))^2,4) 
  } else {
    h2 = h_sq
  }
  
  delta.func = function(x.point){
    value = as.numeric(t(x.point)%*%solve(t(X)%*%X)%*%x.point)
    return(value)
  }
  
  k_factor = rep(NA,length(x0[,1]))
  
  for (j in 1:length(x0[,1])){
    
    x.point=x0[j,]
    delta = delta.func(x.point)
    
    g = (r0^2) - h2 - delta*(c-2*v*k^2)
    
    A = 4*(r0^2) - 4*g
    B = 16*r0*v*k*delta
    C = -8*v*delta*g + 16*(v^2)*(k^2)*(delta^2)
    
    if (B^2 - 4*A*C < 0){　# If the inside of the sqrt is negative
      warning("B^2-$AC is less than 0. It is treated as 0")
      k_factor[j] = round(max(c((-B)/(2*A),(-B)/(2*A))),2)
    } else if (B^2 - 4*A*C >= 0){
      k_factor[j] = round(max(c((-B - sqrt(B^2 - 4*A*C))/(2*A),(-B + sqrt(B^2 - 4*A*C))/(2*A))),2)
    }
    
  } # end for loop
  
  k_factor = abs(k_factor)
  
  output = list(k_factor,h2,c)
  names(output) = c("kfactors","h2","cm")
  return(output)
  
}

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
  
  # mmm = rep(NA,1000)
  # 
  # for (i in 1:1000){
  #   mmm[i]=M_value(X=X,alpha=alpha)  
  #   print(i)
  # } # end for (i in 1:1000)
  # 
  # m = mean(mmm)
  
  
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


# MW_k_factor(P.target=0.75, 
#             alpha.target=0.01,
#             X=X, 
#             x0=x0, 
#             cm=NULL, h_sq=NULL)

###############################################################
# Function to find cm used in modified Wilson's procedure.

# find_cm =function(gamma,X,min_cm,max_cm,cm_by){
#   
#   width =0.0001
#   n = dim(X)[1]
#   q = dim(X)[2]
#   v = n-q
#   k = sqrt((2*v-1)/(2*v))
#   infinity = v+10*sqrt(v)
#   
#   area_under_curve = function(cm,width,X){
#     
#     width = width
#     cm = cm
#     w1 = seq(v*k^2,infinity,by=width)
#     w2 = seq(v*(k-sqrt(cm/(2*v)))^2,v*k^2,by=width)
#     
#     G1 = function(w1) pchisq(q=w1*cm/(v*k^2),df=q)*dchisq(x=w1,df=v)
#     G2 = function(w2) pchisq(q=cm-2*v*(sqrt(w2/v)-k)^2,df=q)*dchisq(x=w2,df=v)
#     
#     return(sum(G1(w1)*width) + sum(G2(w2)*width))
#   } # end of the function area_under_curve
#   
#   Vectorized_area_under_curve = Vectorize(area_under_curve,　vectorize.args = "cm")
#   
#   min_cm = min_cm
#   max_cm = max_cm
#   
#   # create a vector of cm candidate values
#   cm.cand = seq(min_cm,max_cm,by=cm_by)
#   # Find the answer to (2.24) from Limam's dissertation for each cm candidate
#   area = Vectorized_area_under_curve(cm=cm.cand,width = width, X=X )
#   # distance between gamma and (2.24)
#   distance = abs(area - gamma)
#   # plot the distances against cm candidate values
#   plot(cm.cand,distance,type="l")
#   # cm is the value that minimizes the distance.
#   cm = cm.cand[which.min(distance)]
#   
#   return(cm)
# }

# Find cm using integrate and uniroot function
# This is much more efficient!

# find_cm2 =function(gamma,X){
#   
#   n = dim(X)[1]
#   q = dim(X)[2]
#   v = n-q
#   k = sqrt((2*v-1)/(2*v))
#  
#   difference = function(c){
#     cm=c
#     func1 = function(w1) pchisq(q=w1*cm/(v*k^2),df=q)*dchisq(x=w1,df=v)
#     func2 = function(w2) pchisq(q=cm-2*v*(sqrt(w2/v)-k)^2,df=q)*dchisq(x=w2,df=v)
#     G1 = integrate(func1,lower = v*k^2, upper = Inf)$value
#     G2 = integrate(func2,lower = v*(k-sqrt(cm/(2*v)))^2, upper = v*k^2)$value
#   return(G1+G2-gamma)
#   }
#   
#   differenceV = Vectorize(difference)
#   
#   curve(differenceV,xlim=c(0,100))
#   abline(h=0,col="red")
#   
#   return(round(uniroot(differenceV,interval = c(0,50))$root,3))
# 
# }
# 
# find_cm2(gamma=0.90,X=X)
# #[1] 5.018
# find_cm2(gamma=0.95,X=X)
# #[1] 6.432
# find_cm2(gamma=0.99,X=X)
# #[1] 9.653

# D = data.frame(
#   x = c(1.31,1.313,1.32,1.322,1.338,1.34,1.347,1.355,1.36,1.364,1.373,1.376,1.384,1.395,1.4),
#   y = c(4360,4590,4520,4770,4760,5070,5230,5080,5550,5390,5670,5490,5810,6060,5940)
# )
# 
# # Print the observed data table as latex format
# mydata = cbind(rep(1,length(D$x)),D)
# colnames(mydata) = c("intercept","X","Y")
# #print(xtable(mydata,digits=c(0,0,3,0)),include.rownames=FALSE)
# 
# cbind(X0 = rep(1,dim(D)[1]),D)
# 
# X = cbind(rep(1,length(D$x)),D$x)
# x.bar = 1.3531
# N = dim(D)[1]
# x.values = c(1.31,1.353,1.40)
# Sx = sd(D$x)
# 
# z = c(0,0.5,1.0,1.474,1.5,1.604,2.0,2.5,3.0)
# 
# x1 = c(1, 0*Sx + x.bar); x2 = c(1, 0.5*Sx + x.bar) ; x3 = c(1, 1*Sx + x.bar)
# x4 = c(1, 1.474*Sx + x.bar); x5 = c(1, 1.5*Sx + x.bar); x6 = c(1, 1.604*Sx + x.bar)
# x7 = c(1, 2.0*Sx + x.bar); x8 = c(1, 2.5*Sx + x.bar); x9 = c(1, 3.0*Sx + x.bar)
# 
# x0 = rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)
# 
# MW_k_factor(P.target=0.9, alpha.target=0.01, X = X, x0 = x0, cm=NULL, h_sq=NULL)
# # $cm [1] 9.653
# 
# MW_k_factor(P.target=0.9, alpha.target=0.05, X = X, x0 = x0, cm=NULL, h_sq=NULL)
# # $cm [1] 6.432
# 
# xx = x0[c(1,5,9),]
# MW_k_factor(P.target=0.9, alpha.target=0.05, X = xx, x0 = xx, cm=NULL, h_sq=NULL)
