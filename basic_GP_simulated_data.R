# generate true failure data
set.seed(48)
n = 150
x_1 = rnorm(n,2,1) #E = 2
x_2 = rexp(n,1) #E = 1
x_3 = rbinom(n, 1, 0.3) #E = 0.3
x_4 = rnorm(n,-4,7) #E = -4
x_5 = rgamma(n,3,3) #E = 1
x_6 = rpois(n,2) #E = 2
x_7 = asin(runif(n,-1,1)) # E = 0
x_8 = rbeta(n,2,1.5) # E = 0.57
error = rnorm(n,0,0.2)
z = -(x_1^3)/40-(x_1*x_8)/2-x_8-x_2 +1.5
t = exp(error+z)
hist(t,breaks = 50)
# generate the observation time and type of censor
time_length = max(t)
obs_num = 1+rpois(n,10)
obs_time = matrix(0,n,3)
for(i in 1:n)
  {
  obs = rexp(obs_num[i],7)
  left = 0
  j = 1
  right = obs[j]
  type = 1
  while(right < t[i]){
    type = 2
    if(j < obs_num[i]){ 
      left = right
      right = left + obs[j+1]
      j = j+1
    }
    else{
      left = right
      right = Inf
      type = 3
    }  
  }
  obs_time[i,1] = left
  obs_time[i,2] = right
  obs_time[i,3] = type
}


obs_time = as.data.frame(obs_time)
colnames(obs_time) = c('left','right','delta')
obs_time = as.matrix(cbind(obs_time, x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8))
hist(obs_time[,"left"],breaks=50)
hist(obs_time[,'right'],breaks=50)
write.table(obs_time, file="obs_time.csv", sep = " ")
