#total = rpkms$total_introns
#for(i in 1:100){
  #plot(c(0, 0.17, 0.33, 0.5, 1, 2, 4, 8, 16),total[i,])
#  print(i)
  #Sys.sleep(1)
  
#}

library(dlm)
library(ggplot2)



num_data = 1000
t_time = seq(0.01,10,0.01)

# 111, a=3t,b=2,c=1
P_data_true = exp(-t_time)+3*(t_time - 1)
T_data_true = exp(-2*t_time)+3*(6*t_time-7)/4+2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - P_data_true + rnorm(num_data,0,1)
t_dot = 3*t_time - 2*(T_data_true-P_data_true) + rnorm(num_data,0,1)

# 222 a=3,b=2,c=1
P_data_true = exp(-t_time)+3
T_data_true = exp(-2*t_time)+9/2 + 2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3 - P_data_true + rnorm(num_data,0,1)
t_dot = 3 - 2*(T_data_true-P_data_true) + rnorm(num_data,0,1)

# 333 a=3t,c=2,b=1
P_data_true = exp(-2*t_time)+3/4*(2*t_time-1)
T_data_true = exp(-t_time)+3/4*(6*t_time-7)-1
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - 2*P_data_true + rnorm(num_data,0,1)
t_dot = 3*t_time - (T_data_true-P_data_true) + rnorm(num_data,0,1)


truedata=matrix(c(P_data_true,T_data_true), nrow = num_data, byrow=FALSE)
Z = matrix(c(p_dot,t_dot), nrow = num_data, byrow=FALSE)

GG <- diag(1,3)
FF <- matrix(c(1,4,3,4,2,3),nrow=2,byrow=TRUE)
#FF <- matrix(c(2,5,4,2,3,5),nrow=2,byrow=TRUE)
m0 <- c(0,0,0)
C0 <- diag(1,3)
W <- diag(1,3)
V <- diag(1,2)
X = matrix(c(-P_data,P_data-T_data,rep(1,num_data),rep(0,num_data)),nrow=num_data,byrow=FALSE)
#X = matrix(c(T_data,-P_data,-T_data,rep(1,num_data),rep(0,num_data)),nrow=num_data,byrow=FALSE)
my_dlm <- dlm(X=X,JFF=FF,FF=FF,V=V,GG=GG,W=W,m0=m0, C0=C0)


y <- dlmSmooth(Z,my_dlm)
pre_Z = matrix(c(0,0),nrow = 1)
for(k in 1:num_data){
  pre_Z=rbind(pre_Z,t(matrix(c(-P_data[k],0,1,0,(P_data[k]-T_data[k]),1),nrow=2,byrow=TRUE)%*%as.matrix(y$s[k+1,],nrow=3)))
  #pre_Z=rbind(pre_Z,t(matrix(c(-P_data[k],0,1,P_data[k],-T_data[k],0),nrow=2,byrow=TRUE)%*%as.matrix(y$s[k+1,],nrow=3)))
}
pre_Z = dropFirst(pre_Z)

df <- data.frame(y=Z, s=dropFirst(y$s),pre=pre_Z, truedata=truedata)
df$id <- t_time
ggplot(data=df,aes(x=id,y=truedata.1))+ # P
  geom_line()+
  geom_line(aes(y=truedata.2),color="red") # T
  

ggplot(data=df,aes(x=id,y=y.1))+ # P
  geom_point()+
  geom_point(aes(y=y.2),color="red")+ # T
  geom_line(linetype=1,aes(y=pre.1),color="orange")+ #pre_P
  geom_line(linetype=1,aes(y=pre.2),color="darkred") # pre_T

ggplot(data=df,aes(x=id,y=s.1))+ # P
  geom_line()+
  geom_line(linetype=1,aes(y=s.2),color="darkred")+ # b
  geom_line(linetype=1,aes(y=s.3),color="green") #a


num_data = 100
t_time = seq(0.1,10,0.1)

# 111, a=3t,b=2,c=1
P_data_true = exp(-t_time)+3*(t_time - 1)
T_data_true = exp(-2*t_time)+3*(6*t_time-7)/4+2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - P_data_true + rnorm(num_data,0,1)
t_dot = 3*t_time - 2*(T_data_true-P_data_true) + rnorm(num_data,0,1)

# 222 a=3,b=2,c=1
P_data_true = exp(-t_time)+3
T_data_true = exp(-2*t_time)+9/2 + 2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,0.1)
T_data =T_data_true + rnorm(num_data,0,0.1)
p_dot = 3 - P_data_true + rnorm(num_data,0,0.1)
t_dot = 3 - 2*(T_data_true-P_data_true) + rnorm(num_data,0,0.1)

# 333 a=3t,c=2,b=1
P_data_true = exp(-2*t_time)+3/4*(2*t_time-1)
T_data_true = exp(-t_time)+3/4*(6*t_time-7)-1
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - 2*P_data_true + rnorm(num_data,0,1)
t_dot = 3*t_time - (T_data_true-P_data_true) + rnorm(num_data,0,1)


# 444 a=3,c=1/t,b=1/t
P_data_true = 1/t_time +3*t_time/2
T_data_true = 1/t_time + 3*t_time/4 + log(t_time)/t_time
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3 - P_data_true/t_time + rnorm(num_data,0,1)
t_dot = 3 - (T_data_true-P_data_true)/t_time + rnorm(num_data,0,1)



true_a = 3 #rep(3,num_data) #
truedata=matrix(c(P_data_true,T_data_true), nrow = num_data, byrow=FALSE)
Z = matrix(c(p_dot-true_a, t_dot-true_a), nrow = num_data, byrow=FALSE)

GG <- diag(1,2)
FF <- matrix(c(1, 3, 3, 2),nrow=2,byrow=TRUE)
m0 <- c(0,0)
C0 <- diag(1,2)
W <- diag(1,2)
V <- diag(1,2)
X = matrix(c(-P_data,P_data-T_data, rep(0,num_data)),nrow=num_data,byrow=FALSE)
my_dlm <- dlm(X=X,JFF=FF,FF=FF,V=V,GG=GG,W=W,m0=m0, C0=C0)


y0 <- dlmFilter(Z,my_dlm)

y <- dlmSmooth(y0)
#y$s[,1] = tsSmooth(StructTS(x=y$s[,1]))[,1]
#y$s[,2] = tsSmooth(StructTS(x=y$s[,2]))[,1]

pre_Z = matrix(c(0,0),nrow = 1)
for(k in 1:num_data){
  pre_Z=rbind(pre_Z,t(matrix(c(-P_data[k],0,0,(P_data[k]-T_data[k])),nrow=2,byrow=TRUE)%*%as.matrix(y$s[k+1,],nrow=2)))
}
pre_Z = dropFirst(pre_Z)



df <- data.frame(y=Z, s=dropFirst(y$s),pre=pre_Z, truedata=truedata)
df$id <- t_time
ggplot(data=df,aes(x=id,y=truedata.1))+ # P
  geom_line()+
  geom_line(aes(y=truedata.2),color="red") # T


ggplot(data=df,aes(x=id,y=y.1))+ # P
  geom_point()+
  geom_point(aes(y=y.2),color="red")+ # T
  geom_line(linetype=1,aes(y=pre.1),color="orange")+ #pre_P
  geom_line(linetype=1,aes(y=pre.2),color="darkred") # pre_T



ggplot(data=df,aes(x=id,y=s.1))+ # P
  geom_line()+
  geom_line(linetype=1,aes(y=s.2),color="darkred")+ # b
  geom_line(linetype=1,aes(y=1/t_time),color="green")

###############################################

num_data = 100
t_time = seq(0.1,10,0.1)

# 111, a=3t,b=2,c=1
P_data_true = exp(-t_time)+3*(t_time - 1)
T_data_true = exp(-2*t_time)+3*(6*t_time-7)/4+2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - P_data_true + rnorm(num_data,0,1)
t_dot = 3*t_time - 2*(T_data_true-P_data_true) + rnorm(num_data,0,1)

# 222 a=3,b=2,c=1
P_data_true = exp(-t_time)+3
T_data_true = exp(-2*t_time)+9/2 + 2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3 - P_data_true + rnorm(num_data,0,1)
t_dot = 3 - 2*(T_data_true-P_data_true) + rnorm(num_data,0,1)

# 333 a=3t,c=2,b=1
P_data_true = exp(-2*t_time)+3/4*(2*t_time-1)
T_data_true = exp(-t_time)+3/4*(6*t_time-7)-1
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - 2*P_data_true #+ rnorm(num_data,0,1)
t_dot = 3*t_time - (T_data_true-P_data_true) #+ rnorm(num_data,0,1)

# 444 a=3,c=1/t,b=1/t
P_data_true = 1/t_time +3*t_time/2
T_data_true = 1/t_time + 3*t_time/4 + log(t_time)/t_time
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3 - P_data_true/t_time #+ rnorm(num_data,0,1)
t_dot = 3 - (T_data_true-P_data_true)/t_time #+ rnorm(num_data,0,1)


true_a = rep(3,num_data) #3*t_time#
truedata=matrix(c(P_data_true,T_data_true), nrow = num_data, byrow=FALSE)
truedot = matrix(c(p_dot,t_dot), nrow = num_data, byrow=FALSE)
Z = matrix(c(true_a,true_a), nrow = num_data, byrow=FALSE)


GG <- matrix(rep(4,16),nrow=4,byrow=TRUE)
GG[3,3] = 5
GG[4,4] = 6
FF <- matrix(c(1,4,3,4,4,2,4,3),nrow=2,byrow=TRUE)
m0 <- c(0,0,0,0)
C0 <- diag(1,4)
W <- matrix(rep(4,16),nrow=4,byrow=TRUE)
W[3,3] = 7
W[4,4] = 8
V <- diag(1,2)
M_data = T_data-P_data
X = matrix(c(P_data,(T_data-P_data),rep(1,num_data),rep(0,num_data),
             P_data/c(P_data[1],P_data[(length(P_data)-1)]),
             M_data/c(M_data[1],M_data[(length(M_data)-1)]),
             true_a - P_data/c(P_data[1],P_data[(length(P_data)-1)]) * c(true_a[1],true_a[(length(true_a)-1)]) - P_data,
             true_a - M_data/c(M_data[1],M_data[(length(M_data)-1)]) * c(true_a[1],true_a[(length(true_a)-1)]) - M_data
             ),nrow=num_data,byrow=FALSE)
my_dlm <- dlm(X=X,JFF=FF,FF=FF,V=V, JGG=GG, GG=GG, JW=W,W=W,m0=m0, C0=C0)


y <- dlmMLE(Z,my_dlm)
y = dlmSmooth(y)
pre_Z = matrix(c(0,0),nrow = 1)
for(k in 1:num_data){
  pre_Z=rbind(pre_Z,t(matrix(c(P_data[k],0,1,0,0,(T_data[k]-P_data[k]),0,1),nrow=2,byrow=TRUE)%*%as.matrix(y$s[k+1,],nrow=4)))
}
pre_Z = dropFirst(pre_Z)

df <- data.frame(y=Z, s=dropFirst(y$s),pre=pre_Z, truedata=truedata, truedot=truedot)
df$id <- t_time
ggplot(data=df,aes(x=id,y=truedata.1))+ # P
  geom_line()+
  geom_line(aes(y=truedata.2),color="red") # T


ggplot(data=df,aes(x=id,y=y.1))+ # P
  geom_point()+
  geom_point(aes(y=y.2),color="red")+ # T
  geom_line(linetype=1,aes(y=pre.1),color="orange")+ #pre_P
  geom_line(linetype=1,aes(y=pre.2),color="darkred") # pre_T

ggplot(data=df,aes(x=id,y=s.1))+ # c
  geom_line()+
  geom_line(linetype=1,aes(y=s.2),color="darkred")+ # b
  geom_line(linetype=1,aes(y=s.3),color="green")+ #p_dot_pre
  geom_line(linetype=1,aes(y=s.4),color="blue")+ #t_dot_pre
  geom_line(linetype=1,aes(y=truedot.1),color="brown")+ #t_dot
  geom_line(linetype=1,aes(y=truedot.2),color="yellow") #t_dot






###############################################

num_data = 100
t_time = seq(0.1,10,0.1)

# 111, a=3t,b=2,c=1
P_data_true = exp(-t_time)+3*(t_time - 1)
T_data_true = exp(-2*t_time)+3*(6*t_time-7)/4+2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - P_data_true + rnorm(num_data,0,1)
t_dot = 3*t_time - 2*(T_data_true-P_data_true) + rnorm(num_data,0,1)

# 222 a=3,b=2,c=1
P_data_true = exp(-t_time)+3
T_data_true = exp(-2*t_time)+9/2 + 2*exp(-t_time)
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3 - P_data_true + rnorm(num_data,0,1)
t_dot = 3 - 2*(T_data_true-P_data_true) + rnorm(num_data,0,1)

# 333 a=3t,c=2,b=1
P_data_true = exp(-2*t_time)+3/4*(2*t_time-1)
T_data_true = exp(-t_time)+3/4*(6*t_time-7)-1
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3*t_time - 2*P_data_true #+ rnorm(num_data,0,1)
t_dot = 3*t_time - (T_data_true-P_data_true) #+ rnorm(num_data,0,1)

# 444 a=3,c=1/t,b=1/t
P_data_true = 1/t_time +3*t_time/2
T_data_true = 1/t_time + 3*t_time/4 + log(t_time)/t_time
P_data = P_data_true + rnorm(num_data,0,1)
T_data =T_data_true + rnorm(num_data,0,1)
p_dot = 3 - P_data_true/t_time #+ rnorm(num_data,0,1)
t_dot = 3 - (T_data_true-P_data_true)/t_time #+ rnorm(num_data,0,1)


true_a = rep(3,num_data) #3*t_time#

truedata=matrix(c(P_data_true,T_data_true), nrow = num_data, byrow=FALSE)
truedot = matrix(c(p_dot,t_dot), nrow = num_data, byrow=FALSE)

true_a_t_1 = c(true_a[1],true_a[1:(length(true_a)-1)])
P_data_t_1 = c(P_data[1],P_data[1:(length(P_data)-1)])
Z = matrix(true_a, nrow = num_data, byrow=FALSE)


GG <- matrix(c(3,4,10,3),nrow=2,byrow=TRUE)
FF <- matrix(c(1,3),nrow=1,byrow=TRUE)
m0 <- c(0,0)
C0 <- diag(1,2)
w = 0.1
W <- matrix(c(7,8,8,9),nrow=2,byrow=TRUE)
V <- diag(1,1)
M_data = T_data-P_data
X = matrix(c(P_data,(T_data-P_data),rep(1,num_data),rep(0,num_data),
             P_data/c(P_data[1],P_data[(length(P_data)-1)]),
             M_data/c(M_data[1],M_data[(length(M_data)-1)]),
             rep(w**2, num_data),
             w*(true_a - true_a_t_1 + w*P_data),
             (true_a - true_a_t_1 + w*P_data)**2,
             P_data_t_1 - P_data
),nrow=num_data,byrow=FALSE)

my_dlm <- dlm(X=X,JFF=FF,FF=FF,V=V, JGG=GG, GG=GG, JW=W,W=W,m0=m0, C0=C0)


y <- dlmFilter(Z,my_dlm)
y = dlmSmooth(y)
pre_Z = matrix(c(0),nrow = 1)
for(k in 1:num_data){
  pre_Z=rbind(pre_Z,t(matrix(c(P_data[k],1),nrow=1,byrow=TRUE)%*%as.matrix(y$s[k+1,],nrow=2)))
}
pre_Z = dropFirst(pre_Z)

df <- data.frame(y=Z, s=dropFirst(y$s),pre=pre_Z, truedata=truedata, truedot=truedot)
df$id <- t_time
ggplot(data=df,aes(x=id,y=truedata.1))+ # P
  geom_line()+
  geom_line(aes(y=truedata.2),color="red") # T


ggplot(data=df,aes(x=id,y=y))+ # P
  geom_point()+
  geom_line(linetype=1,aes(y=pre),color="orange") #pre_P

ggplot(data=df,aes(x=id,y=s.1))+ # c
  geom_line()+
  geom_line(linetype=1,aes(y=s.2),color="darkred")+ # b
  geom_line(linetype=1,aes(y=truedot.1),color="brown") # p_dot

