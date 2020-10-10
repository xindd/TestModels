## <--------------------------------------------------------------------------->
## Example 1: ARMA(2, 1) model estimation.
## <--------------------------------------------------------------------------->
## This example shows how to fit an ARMA(2, 1) model using this Kalman
## filter implementation (see also stats' makeARIMA and KalmanRun).
library(RUnit) # 依赖
library(FKF)
num_data = 100
t_time = seq(0.1,10,0.1)

# 444 a=3,c=1/t,b=1/t
P_data_true = 1/t_time +3*t_time/2
T_data_true = 1/t_time + 3*t_time/4 + log(t_time)/t_time
P_data = P_data_true + rnorm(num_data,0,.1)
T_data =T_data_true + rnorm(num_data,0,.1)
p_dot = 3 - P_data_true/t_time #+ rnorm(num_data,0,1)
t_dot = 3 - (T_data_true-P_data_true)/t_time #+ rnorm(num_data,0,1)

true_a = rep(3,num_data)+ rnorm(num_data,0,.1) #3*t_time + rnorm(num_data,0,1)# 
truedata=matrix(c(P_data_true,T_data_true), nrow = num_data, byrow=FALSE)
truedot = matrix(c(p_dot,t_dot), nrow = num_data, byrow=FALSE)

true_a_t_1 = c(true_a[1],true_a[1:(length(true_a)-1)])
P_data_t_1 = c(P_data[1],P_data[1:(length(P_data)-1)])

y = matrix(true_a, nrow = 1, byrow=FALSE)
# dt <- rbind(rep(0, num_data), true_a-true_a_t_1*P_data/P_data_t_1)
dt <- rbind(rep(0, num_data), true_a-true_a_t_1)
ct <- matrix(0)
Zt <- array(rbind(P_data, rep(1, num_data)),dim=c(1,2,num_data))

Tt <- cbind(matrix(rep(1,num_data), nrow=num_data, byrow=TRUE), P_data_t_1 - P_data)
Tt <- cbind(Tt, matrix(rep(0,num_data), nrow=num_data, byrow=TRUE))
Tt <- t(cbind(Tt, matrix(rep(1,num_data), nrow=num_data, byrow=TRUE)))
dim(Tt) <- c(1,4*num_data)
Tt <- array(Tt,dim=c(2,2,num_data))

a0 <- c(0,0)
P0 <- diag(.1,2) # Variance of 'a0'

w <- 0.1
HHt <- cbind(matrix(rep(1,num_data), nrow=num_data, byrow=TRUE), P_data)
HHt <- cbind(HHt, P_data)
HHt <- t(cbind(HHt, P_data**2))
dim(HHt) <- c(1,4*num_data)
HHt <- w**2 * array(HHt,dim=c(2,2,num_data))

v=0.1
GGt = v**2 * array(rep(1,num_data),dim=c(1,1,num_data))

## Filter Nile data with estimated parameters:
fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = HHt,GGt = GGt, yt = y)
plot(fkf.obj, type = "state")
## Plot the flow data :
pre_Z = matrix(c(0),nrow = 1)
for(k in 1:num_data){
  pre_Z=rbind(pre_Z,t(matrix(c(P_data[k],1),nrow=1,byrow=TRUE)%*%as.matrix(fkf.obj$att[,k],nrow=2)))
}
pre_Z = dropFirst(pre_Z)

df <- data.frame(s=t(fkf.obj$att),truedot=truedot, pre_Z = pre_Z, true_a = true_a)
df$id <- t_time
ggplot(data=df,aes(x=id,y=true_a))+ # P
  geom_point()+
  geom_line(linetype=1,aes(y=pre_Z),color="orange") #pre_P

ggplot(data=df,aes(x=id,y=s.1))+ # c
  geom_line()+
  geom_line(linetype=1,aes(y=s.2),color="darkred")+ # pre_dot
  geom_line(linetype=1,aes(y=truedot.1),color="green") # p_dot
