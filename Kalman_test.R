kl <- function(P_data, T_data, p_dot, t_dot, PL_data, t_time, tL){
  #t_time = c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
  #tL = 1/6
  
  num_data = length(t_time)
  true_a = PL_data/tL
  
  
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
  pre_Z = matrix(c(0,0),nrow = 1)
  for(k in 1:num_data){
    pre_Z=rbind(pre_Z,t(matrix(c(-P_data[k],0,0,(P_data[k]-T_data[k])),nrow=2,byrow=TRUE)%*%as.matrix(y$s[k+1,],nrow=2)))
  }
  pre_Z = dropFirst(pre_Z)
  
  df <- data.frame(y=Z, s=dropFirst(y$s),pre=pre_Z)
  df$id <- t_time
  
  ggplot(data=df,aes(x=id,y=y.1))+ # P
    geom_point()+
    geom_point(aes(y=y.2),color="red")+ # T
    geom_line(linetype=1,aes(y=pre.1),color="orange")+ #pre_P
    geom_line(linetype=1,aes(y=pre.2),color="darkred") # pre_T
  
  ggplot(data=df,aes(x=id,y=s.1))+ # P
    geom_line()+
    geom_line(linetype=1,aes(y=s.2),color="darkred") # b
  
}

num_data = 100
t_time = seq(0.1,10,0.1)

P_data = rpkms$total_introns['77595',]
T_data = rpkms$total_exons['77595',]
p_dot = 
t_dot = 
PL_data = rpkms$foursu_exons['77595',]
t_time = c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
tL = 1/6

kl(P_data, T_data, p_dot, t_dot, PL_data, t_time, tL)
