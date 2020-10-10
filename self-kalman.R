
buildFun <- function(Z,num_data,Ft){
  X = matrix(c(0,0), nrow=2)
  P = diag(1,2)
  Q = diag(0.00001,2)
  H = matrix(c(1,0), nrow=1)
  R = 1
  # FF <- matrix(c(1,1,0,1),nrow=2,byrow=TRUE)
  
  FFt <- cbind(matrix(rep(1,num_data), nrow=num_data, byrow=TRUE), rep(0,num_data))
  FFt <- cbind(FFt, Ft)
  FFt <- t(cbind(FFt, rep(1,num_data)))
  dim(FFt) <- c(1,4*num_data)
  FFt <- array(FFt,dim=c(2,2,num_data))
  
  x1 = c()
  x2 = c()
  for(i in 1:num_data){
    FF = FFt[,,1]
    X_x = FF%*%X
    P_p = FF%*%P%*%t(FF) + Q
    K = P_p%*%t(H) / as.numeric((H%*%P_p%*%t(H) + R))
    X = X_x + K%*%(Z[i] - H%*%X_x)
    P = (diag(1,2) - K%*%H)%*%P_p
    
    #FF = matrix(c(1,1+i*0.1,0,1),nrow=2,byrow=TRUE)
    
    x1 = c(x1, X[1])
    x2 = c(x2, X[2])
    
  }
  return(cbind(x1,x2))
}


t_time = c(0, 0.17, 0.33, 0.5, 1, 2, 4, 8, 16)
dert = c(1/6, diff(t_time))
#tessdata = c(2.852587, 2.477313, 2.732977, 2.333752, 2.117182, 2.749777, 2.334962, 2.472333, 2.129490)
tessdata = c(2.852587, 2.477313, 2.732977, 2.333752, 2.117182, 2.749777, 2.334962, 2.472333, 2.129490)
num_data = length(t_time)

Z = tessdata
res = buildFun(Z, num_data, Ft=dert)

df <- data.frame(s=res, truedata=tessdata)
df$id <- t_time

ggplot(data=df,aes(x=id,y=s.x1))+ # pre_data
  geom_line(color="black")+
  geom_line(linetype=1,aes(y=s.x2),color="darkred")+ # pre_dot
  geom_point(aes(y=truedata),color="darkred") # pre_data


     