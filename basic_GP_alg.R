# check packages
check_packages = function(names)
{
  for(name in names)
  {
    if (!(name %in% installed.packages()))
      install.packages(name, repos="http://cran.us.r-project.org")
    
    library(name, character.only=TRUE)
  }
}

check_packages(c("PCDSpline","utils","truncnorm","MASS","MCMCpack",'coda','devtools'))

#Data import
data = read.csv("obs_time.csv", header = T, sep = " ") 
data_row = nrow(data) # number of observations
x.train = as.matrix(data[,4:ncol(data)])
var_col = ncol(data)-3 # number of predictors 
delta = as.data.frame(model.matrix(~ factor(data$delta) - 1)) #factorize indicators for censoring
colnames(delta) = c('delta 1','delta 2','delta 3')
delta_1 = as.data.frame(model.matrix(~ factor(delta$`delta 1`) - 1))
t = data[c('left','right')]*delta_1  
t = rowSums(t, na.rm = TRUE) #matrix of t
delta_2_index = which(delta['delta 2']==1) #row index whose delta 2 is 1

# number of basis functions and generating knots
k = 20  #20+1 knots
degree = 4 
time_span = data[c('left','right')]
time_span = do.call(data.frame,lapply(time_span, function(x) replace(x, is.infinite(x),NA)))
time_span = do.call(data.frame,lapply(time_span, function(x) replace(x, x == 0,NA)))
ksi = seq(from = 0, to = max(time_span,na.rm = TRUE)/k*(k+1),length.out = k+1) #does the knot placement matters?

# generating the basis functions
ispline_t = Ispline(t, degree, ksi)
ispline_L = Ispline(data[,'left'],degree, ksi) #(alpha(L))
ispline_R = Ispline(data[,'right'],degree, ksi) #(alpha(R))

#prior parameters

v0 = 1/4 #prior precision of alpha_0
m0 = 0 #prior mean of alpha_0
a = 2 #prior shape of eta
b = 2 #prior rate of eta
v_a = 4 #prior shape of v, the variance of signal
v_b = 5 #prior rate of v, the variance of signal
a_1 = 2 #prior shape of eta_d
b_1 = 2 #prior rate of eta_d
a_2 = 2 #prior shape of b_d
b_2 = 2 #prior rate of b_d
l = sqrt(var_col)
#l = 1

# Initialization
iter = 3000 #iteration times

alpha_0 = matrix(0,nrow = iter + 1) #save the constant term in the spline regression = gamma_0 in write-up
alpha = matrix(0.1, ncol = k+degree-1, nrow = iter + 1) #save the coefficient = gamma_I in write-up
f = matrix(0,ncol = data_row,nrow = iter + 1) #save values of f(x), in each iteration
z = matrix(1, ncol = data_row, nrow=1) #temporary vector that contains z (data argumentation)
eta = matrix(1,nrow = iter+1) #hyperparamter of spline coefficients 
v = matrix(1, nrow = iter+1) #variance of signal
inv_length_scale = matrix(1, ncol = var_col, nrow = iter + 1) #inverse length scale, saved in each iteration
inv_length_scale_temp = diag(inv_length_scale[1,]) #temporary diagonal matrix of inverse length scale, used in each iteration
eta_d = matrix(1, nrow = iter+1) #hyperparameter of length-scale
b_d = matrix(1, nrow = iter+1) #hyperparameter of length-scale
#eta_d = matrix(1, nrow = iter+1)

# function to calculate S
S_func = function(X, l_s,l){ #l revise
  #each row of X is an observation
  #l_s must be a diagonal matrix
  #n is number of observations
  #can add a step to check the input
  n = nrow(X)
  S = diag(0,n)
  for(a in 1:n){
    for(b in a:n){  #symmetric covariance matrix
      S[a,b]=S[b,a]=exp(-1/2*(X[a,]-X[b,])%*%(l_s^l)%*%t(t(X[a,]-X[b,])))
    }
  }
  return(S)
}

# calculate the initial S (S_old)
S_old = S_func(x.train, inv_length_scale_temp,l)
inv_S_old = solve(S_old)
  
#MCMC
for (j in 1:iter){
  cat(j,",")
  #sampling z
  for(i in 1:data_row){
    z_mean = alpha_0[j] + alpha[j,]%*%ispline_t[,i] + f[j,i] #check sum or colsum
    if(delta[i,'delta 1']==1){
      z[i] = rtruncnorm(1, a = 0, b = Inf, mean = z_mean, sd = 1)
    } else if(delta[i,'delta 2']==1){
      low = alpha[j,]%*%ispline_L[,i] - alpha[j,]%*%ispline_R[,i] #check sum or colsum
      z[i] = rtruncnorm(1, a = low, b = 0, mean = z_mean, sd = 1)
    } else{
      z[i] = rtruncnorm(1, a=-Inf, b = 0, mean = z_mean, sd = 1)
    }
  }
  
  #sample alpha_0
  W_0 = v0 + data_row
  E_0 = (v0*m0 + sum(z - alpha[j,]%*%ispline_t - f[j,]))/W_0 #check sum or colsum
  alpha_0[j+1]= rnorm(1,E_0,(1/W_0)^0.5)
  
  #sample alpha_I
  for(I in 1:(k+degree-1)){
    W = sum(ispline_t[I,]^2)
    if(W==0){
      alpha[j+1,I]=rexp(1,eta[j])
    }
    else{
      c_I = max((-z[delta_2_index]-alpha[j,-I]%*%
                   (ispline_R[-I,delta_2_index]-ispline_L[-I,delta_2_index]))
                /(ispline_R[I,delta_2_index]-ispline_L[I,delta_2_index]))
      # not exactly gibbs in the above step, but the MC should be ergodic 
      d_I = max(c_I,0)
      E_I = ((z-alpha_0[j+1]-alpha[j,-I]%*%ispline_t[-I,]-f[j,])%*%t(t(ispline_t[I,]))-eta[j])/W
      alpha[j+1,I]=rtruncnorm(1, a = d_I, b = Inf, mean = E_I, sd = (1/W)^0.5)
    }
  } #recheck required
  
  #sampling eta
  eta[j+1] = rgamma(1, a+k+degree-1, b+sum(alpha[j+1,]))
  
  #sampling f
  spline = alpha_0[j+1] + alpha[j+1,]%*%ispline_t #subset may change horizontal/vertical
  W_f = diag(data_row) + solve(v[j]*S_old) #100*100
  iW_f = solve(W_f)
  E_f = iW_f%*%t(z-spline) #might have risk, check vertical or horizontal 100*1
  f[j+1,] = mvrnorm(1,E_f,iW_f)
  
  #sampling v
  v[j+1] = rinvgamma(1, v_a + data_row/2, v_b + 0.5*f[j+1,]%*%inv_S_old%*%t(t(f[j+1,])))

  #sampling m, the inverse length-scale
  inv_length_scale_temp = diag(inv_length_scale[j,]) #extract length-scale values in the last iteration and make a diagonal matrix
  for(p in 1:var_col){
    prop = rtruncnorm(1, a=0, b=Inf, mean = inv_length_scale[j,p], sd = 1)
    inv_length_scale_temp[p,p] = prop

    #calculate S
    S = S_func(x.train, inv_length_scale_temp,l)
    inv_S = solve(S)
    
    #calculate acceptance ratio of each length-scale
    alpha_num = -0.5*f[j+1,]%*%inv_S%*%t(t(f[j+1,]))/v[j+1]
    - 1/2*log(det(S*v[j+1])) + (eta_d[j]*l-1)*log(prop) - 
      b_d[j]*prop^l - log(dtruncnorm(prop, a=0, b=Inf, mean = inv_length_scale[j,p], sd = 1))
    
    alpha_den = -0.5*f[j+1,]%*%inv_S_old%*%t(t(f[j+1,]))/v[j+1]
    + 1/2*log(det(inv_S_old/v[j+1])) + (eta_d[j]*l-1)*log(inv_length_scale[j,p])- 
      b_d[j]*inv_length_scale[j,p]^l-log(dtruncnorm(inv_length_scale[j,p], a=0, b=Inf, mean = prop, sd = 1))
    accept = min(1, exp(alpha_num-alpha_den)) 

    U = runif(1,0,1)
    if(U > accept){
      inv_length_scale_temp[p,p] = inv_length_scale[j,p] #reject
      #S_old = S_func(x.train, inv_length_scale_temp,l)
      #inv_S_old = solve(S_old)
    }else{
      S_old = S
      inv_S_old = inv_S
    }
    #cat(inv_length_scale_temp)
  }
  inv_length_scale[j+1,] = diag(inv_length_scale_temp)
  
  S_old = S_func(x.train, inv_length_scale_temp,l) #redundant
  inv_S_old = inv_S #redundant
  
  #sample eta_d, M_H alg
  prop = rtruncnorm(1, a=0, b=Inf, mean = eta_d[j], sd = 1)
  alpha_num = (a_1-1)*log(prop)-b_1*prop+var_col*(prop*log(b_d[j])-log(gamma(prop)))
  +(prop-1)*l*sum(log(inv_length_scale[j+1,]))-log(dtruncnorm(prop, a=0, b=Inf, mean = eta_d[j], sd = 1))
  
  alpha_den = (a_1-1)*log(eta_d[j])-b_1*eta_d[j]+var_col*(eta_d[j]*log(b_d[j])-log(gamma(eta_d[j])))
  +(eta_d[j]-1)*l*sum(log(inv_length_scale[j+1,]))-log(dtruncnorm(eta_d[j], a=0, b=Inf, mean = prop, sd =1))
  
  accept = min(1, exp(alpha_num-alpha_den))
  
  U = runif(1,0,1)
  if(U < accept){
    eta_d[j+1]=prop 
  } else{
    eta_d[j+1]=eta_d[j]
  }
  
  
  #sample b_d
  b_d[j+1] = rgamma(1,var_col*eta_d[j+1]+a_2, b_2+sum(inv_length_scale[j+1,]^l)) #real inv-scale = inv^l
  
  #ncomp: no. of components
  #pmtm.budget = 5*ncomp
  #p.incl: prob of inclusion
  #varp.update: variable importance score vector? why initialized in 99, but not 1?
  #lpentalty: log penalty score for knots selection 
  #out = add.gp(x.train, spline, nsweep = 1e3, prox = c(0.999, 0.99, 0.9, 0.7), Rsq = c(0.2, 0.6, 0.99, 0.999), lp.lam.f = 0, lp.rho.f = 0, ncomp, p.incl, p.ref = 0.1, a = 1, b = 1, max.rank = nrow(x.train), tol = c(1e-4,1e-6), burn = 0.1, pmtm.budget, varp.update, lpenalty)
  #out = add.gp(x.train, t(spline), nsweep = 1e3, prox = c(0.99, 0.94, 0.88, 0.81, 0.73), Rsq = c(0.01, 0.25, 0.5, 0.7, 0.86, 0.99), lp.lam.f = 0, lp.rho.f = 0, ncomp, p.incl, p.ref = 0.1, a = 1, b = 1, max.rank = nrow(x.train), tol = c(1e-4,1e-6), burn = 0.1, pmtm.budget, varp.update, lpenalty)
}


#normality test
burn = 1
alpha_0_b = as.matrix(alpha_0[burn:iter+1,])
plot(alpha_0_b,type= "l")
alpha_b = alpha[burn:iter+1,]
f_b = f[burn:iter+1,]
plot(f_b[,10],type= "l")
plot(alpha_b[,13],type= "l")
inv_length_scale_b = inv_length_scale[burn:iter+1,]
plot(inv_length_scale_b[,1],type= "l")
Z = alpha_0_b %*% rep(1,dim(ispline_t)[2]) + alpha_b %*% ispline_t + f_b #iter*obs
mean_Z = colMeans(Z)
norm_test = as.matrix(shapiro.test(mean_Z))
qqnorm(mean_Z)





write.table(norm_test, file="basic.csv",sep= " ")
write.table(eta_d, file="eta_d.csv", sep= " ")
write.table(f_b, file="f.csv", sep= " ")
write.table(inv_length_scale_b, file="inv_length_scale.csv", sep= " ")
write.table(mean_Z, file="residual.csv", sep= " ")
