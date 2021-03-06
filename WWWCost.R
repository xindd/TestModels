library(ggplot2)


wwwcost <- function(N, tL, PL, TL, PR, TR, Length, TimeGrid){
  Stop1 = 1.0
  Stop2 = 1.0

  PError = 0
  PError2 = 0

  alpha = 10000*matrix(1,2,1)
  ####################
  # 最大公约数
  getGCD <- function(m,n){
    if(!n){
      return(m)
    }else{
      return(getGCD(n,m %% n))
    }
  }
  TimeGrid <- c(0, 2,4,6,8)
  TimeGrid2 = TimeGrid[-which(TimeGrid==0)]
  tmp = getGCD(TimeGrid2[1],TimeGrid2[2])
  for(i in 3:length(TimeGrid2)){
    tmp = getGCD(tmp,TimeGrid2[i])
  }
  foldGCD = 2
  tmp = tmp/foldGCD
  gridnumber = round(diff(TimeGrid)/tmp-1)
  timegrid = c()
  for(i in 1:length(gridnumber)){
    if(gridnumber[i]==0){
      timegrid = c(timegrid, TimeGrid[i])
    }else{
      timegrid = c(timegrid, spline(c(TimeGrid[i],(TimeGrid[i+1]-tmp)),n=gridnumber[i]+1)$y)
    }
  }
  #######################
  timegrid = c(timegrid, TimeGrid[i+1])
  pos_at = 1:length(TimeGrid)+c(0,cumsum(gridnumber))
  pos_at_nan = setdiff(1:length(timegrid),pos_at)
  #######################
  TL = labexon_inspect[65,]
  TT = totexon_inspect[65,]
  PT = totintr_inspect[65,]
  N = length(timegrid)
  #######################
  a = matrix(0,1,N)
  b = matrix(0,1,N)
  c = matrix(0,1,N)
  h = tmp
  Time = timegrid
  tL = 1
  # 根据标记数据，设置时间点 at 的值
  a[pos_at] = TL/tL

  #a = seq(a[1],a[length(a)],(a[length(a)]-a[1])/N)
  # 根据标记数据，设置时间点 ct 的值
  #Coef = c(tL*PL[1]/TL[1],tL*PL[length(PL)]/TL[length(TL)])

  #f = function(x){(1 - exp(-tL*x))/x - Coef}
  #g = function(x){(tL*x*exp(-tL*x) - (1 - exp(- tL*x)))/x^2}
  #Temp =  c(c[1],c[length(c)])
  #for(i in  1:10){
  #  Temp = Temp - f(Temp)/g(Temp)
  #}
  #c = seq(Temp[1], Temp[length(Temp)],(Temp[length(Temp)] - Temp[1])/N)

  # 设置矩阵
  if(length(pos_at_nan)==0){
    IA=0
  }else{
    IA = diag(2,length(pos_at_nan),length(pos_at_nan))# +
    #cbind(matrix(0,N-length(pos_at),1),diag(-1,N-length(pos_at), N-length(pos_at)-1)) +
    #rbind(matrix(0,1,N-length(pos_at)),diag(-1,N-length(pos_at)-1, N-length(pos_at)))
    rownames(IA)=pos_at_nan
    colnames(IA)=pos_at_nan
    for(rows in as.character(pos_at_nan)){
      for(cols in as.character(pos_at_nan))
        if(abs(as.numeric(rows)-as.numeric(cols))==1){
          IA[rows, cols]=-1
        }
    }
  }
  IB = diag(2,N,N) + cbind(matrix(0,N,1),diag(-1,N, N-1)) + rbind(matrix(0,1,N),diag(-1,N-1, N))
  IB[1,1] = 1
  IB[dim(IB)[1],dim(IB)[1]] = 1
  IC = IB

  blkdiag <- function(x,y,z){
    x1 = cbind(x,matrix(0, dim(x)[1], dim(y)[1] + dim(z)[1]))
    y1 = cbind(matrix(0, dim(y)[1], dim(x)[1]), y)
    y1 = cbind(y1,matrix(0, dim(y)[1], dim(z)[1]))
    z1 = cbind(matrix(0,dim(z)[1],dim(x)[1]+dim(y)[1]), z)
    mat = rbind(x1,y1)
    mat = rbind(mat, z1)
    return(mat)
  }
  # X 矩阵
  X = c(a[pos_at_nan],b,c)
  # A 矩阵
  MatrixJ = blkdiag(IA,IB,IC)/(2*h)
  # B 矩阵
  VectorJ = matrix(0,3*N-length(pos_at),1)
  for(i in 1:length(pos_at_nan)){
    if(any((pos_at_nan[i]-1)!=pos_at_nan) & any((pos_at_nan[i]+1)!=pos_at_nan)){
      VectorJ[i] = a[pos_at_nan[i]-1] + a[pos_at_nan[i]+1]
    }else if(any((pos_at_nan[i]-1)==pos_at_nan) & any((pos_at_nan[i]+1)==pos_at_nan)){
      VectorJ[i] = 0
    }else if(any((pos_at_nan[i]-1)!=pos_at_nan) & any((pos_at_nan[i]+1)==pos_at_nan)){
      VectorJ[i] = a[pos_at_nan[i]-1]
    }else if(any((pos_at_nan[i]-1)==pos_at_nan) & any((pos_at_nan[i]+1)==pos_at_nan)){
      VectorJ[i] = a[pos_at_nan[i]+1]
    }
  }
  VectorJ = VectorJ/h
  PRTRX <- function(a,b,c,N,timegrid,TT,PT){
    eta = matrix(0,3*N-length(TimeGrid),N)
    for(i in 1:N){
      if(any(i==pos_at_nan)){
        eta[i-1,i] = 1
      }
    }

    xi =  matrix(0,3*N-length(TimeGrid),N)
    for(i in 1:N){
      xi[2*N-length(TimeGrid)+i,i] = 1
    }

    theta = matrix(0,3*N-length(TimeGrid),N)
    for(i in 1 : N){
      theta[N-length(TimeGrid)+i,i] = 1
    }

    # return P_dot, T_dot
    f1 = function(x,a,b,c){
      t(t(c(a - c * x[1],a - b * (x[2]-x[1]))))
    }
    # return P_X_dot, T_X_dot
    g1 = function(x,eta,xi,theta,b,c,S){
      t(t(c(eta - xi * S[1] - c * x[1:(3*N-length(TimeGrid))],
            eta - (S[2]-S[1])*theta - b*(x[(3*N-length(TimeGrid)+1):length(x)]-x[1:(3*N-length(TimeGrid))]))))
    }

    inital_exp = rbind(PT[1], TT[1])
    S=inital_exp
    K = matrix(0,2,4)

    GRAD = matrix(0,2*(3*N-length(TimeGrid)),1)
    L = matrix(0,2*(3*N-length(TimeGrid)),4)
    # 利用龙格库塔计算预测值P，T
    for(i in 1 : (N - 1)){
      S1 = S[,i]
      K[,1] = f1(S1,a[i],b[i],c[i])
      S2 = S[,i] + h/2 * K[,1]
      K[,2] = f1(S2,a[i],b[i],c[i])
      S3 = S[,i] + h/2 * K[,2]
      K[,3] = f1(S3,a[i],b[i],c[i])
      S4 = S[,i] + h * K[,3]
      K[,4] = f1(S4,a[i],b[i],c[i])
      S = cbind(S, S[,i] + h * ( K[,1] + 2 * K[,2] + 2 * K[,3] + K[,4])/6)

      # 利用龙格库塔计算预测值GP，GT
      L[,1] = g1(GRAD, eta[,i], xi[,i], theta[,i], b[i], c[i], S1)
      L[,2] = g1(GRAD+h/2 * L[,1], eta[,i], xi[,i], theta[,i], b[i], c[i], S2)
      L[,3] = g1(GRAD+h/2 * L[,2], eta[,i], xi[,i], theta[,i], b[i], c[i], S3)
      L[,4] = g1(GRAD+h * L[,3], eta[,i], xi[,i], theta[,i], b[i], c[i], S4)

      GRAD = GRAD + h * ( L[,1] + 2 * L[,2] + 2 * L[,3] + L[,4])/6
    }

    dim(GRAD) = c(3*N-length(TimeGrid),2)

    return(list(S,GRAD))

  }
  for(i in 1 : 50){
    # 梯度项1
    Grad = MatrixJ %*% t(t(X)) + VectorJ

    res = PRTRX(a,b,c,N,timegrid,TT,PT)

    S = res[[1]]
    GRAD = res[[2]]
    Grad = Grad + alpha[1] * sum(S[1,pos_at] - PT) * GRAD[,1] + alpha[2] * sum(S[2,pos_at] - TT) * GRAD[,2]
    #Grad = Grad + alpha[1] * (S[1] - PR[length(PR)]) * GRAD[,1] + alpha[2] * (S[2] - TR[length(TR)]) * GRAD[,2]

    Index = which(Grad>0)

    rho = 0.5 * min(X[Index]/Grad[Index])

    X = X - rho * Grad
    StopRule = min(c(Stop1,Stop2,sqrt(sum(Grad^2))))
    print(Grad^2)
    if(StopRule < 1e-4){
      break
    }

    a[pos_at_nan] = X[1:(N-length(TimeGrid))]
    b = X[(N-length(TimeGrid)+1):(2*N)]
    c = X[(2*N+1):(3*N-length(TimeGrid))]

    Error = t(X) %*% MatrixJ %*% X/2/h + t(X) %*% VectorJ/h
    Error2 = (S[1] - PR[length(PR)])^2 + (S[2] - TR[length(TR)])^2
    Stop1 = abs(Error-PError)
    Stop2 = abs(Error2-PError2)
    PError = Error
    PError2= Error2
  }
  Error = t(X) %*% MatrixJ %*% X/2/h + t(X) %*% VectorJ/h
  Error2 = (S[1] - PR[length(PR)])^2 + (S[2] - TR[length(TR)])^2
  return(list(a,b,c,Time))
}
################################
TL = labexon_inspect[65,]/10
TT = totexon_inspect[65,]
PT = totintr_inspect[65,]
TimeGrid <- c(0, 2, 4, 6,8)
tL = 1
res = wwwcost(N, tL, PL, TL, PR, TR, Length, TimeGrid)

################################

TimeGrid = c(0,2)
Length = 2
PL = c(0.1869, 0.3816)
TL = c(0.4898, 0.6463)
PR = c(0.4387, 0.4456)
TR = c(0.7655, 0.7952)
tL = 1.5
N = 30
res = wwwcost(N, tL, PL, TL, PR, TR, Length, TimeGrid)
fa = res[[1]]
fb = res[[2]]
fc = res[[3]]
ftime = res[[4]]
plot(c(0,ftime,2), fa, ylim=c(0,2), type = 'o')
points(c(0,ftime,2), fb, col='red', type = 'o')
points(c(0,ftime,2), fc, col='blue', type = 'o')

getGCD <- function(m,n){
  if(!n){
    return(m)
  }else{
    return(getGCD(n,m %% n))
  }
}
TimeGrid <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
TimeGrid2 = TimeGrid[-which(TimeGrid==0)]
tmp = getGCD(TimeGrid2[1],TimeGrid2[2])
for(i in 3:length(TimeGrid2)){
  tmp = getGCD(tmp,TimeGrid2[i])
}
gridnumber = round(diff(TimeGrid)/tmp-1)
timegrid = c()
for(i in 1:length(gridnumber)){
  if(gridnumber[i]==0){
    timegrid = c(timegrid, TimeGrid[i])
  }else{
    timegrid = c(timegrid, spline(c(TimeGrid[i],(TimeGrid[i+1]-tmp)),n=gridnumber[i]+1)$y)
  }
}
timegrid = c(timegrid, TimeGrid[i+1])
pos_at = 1:length(TimeGrid)+c(0,cumsum(gridnumber))


TimeGrid <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
TL = as.vector(rpkms$foursu_exons[1,])
TR = as.vector(rpkms$total_exons[1,])
PL = as.vector(rpkms$foursu_introns[1,])
PR = as.vector(rpkms$total_introns[1,])
tL = 1/6
Length = 2
N = 30
i=6
timegrid = TimeGrid[c(i,(i+1))]
TL1 = TL[c(i,(i+1))]
TR1 = TR[c(i,(i+1))]
PL1 = PL[c(i,(i+1))]
PR1 = PR[c(i,(i+1))]
res = wwwcost(N, tL, PL1, TL1, PR1, TR1, Length, timegrid)
fa = res[[1]]
fb = res[[2]]
fc = res[[3]]
ftime = res[[4]]


ad = as.matrix(data.frame(ftime=round(ftime,4),aa=round(fa[2:(length(fa)-1)],4), type=rep(1,length(ftime))))
bd = as.matrix(data.frame(ftime=round(ftime,4),bb=round(fb[2:(length(fb)-1)],4), type=rep(2,length(ftime))))
cd = as.matrix(data.frame(ftime=round(ftime,4),cc=round(fc[2:(length(fc)-1)],4), type=rep(3,length(ftime))))
plotdata = rbind(ad,bd)
plotdata = rbind(plotdata, cd)
plotdata = as.data.frame(plotdata)
plotdata[,3] = factor(plotdata[,3])
ggplot(plotdata, aes(x=ftime, y=aa, colour = type, shape = type)) +
  geom_point()# +facet_wrap( ~ type)





TimeGrid <- c(0, 1/6, 1/3, 1/2, 1, 2, 4, 8, 16)
TL = as.vector(rpkms$foursu_exons[1,])
TR = as.vector(rpkms$total_exons[1,])
PL = as.vector(rpkms$foursu_introns[1,])
PR = as.vector(rpkms$total_introns[1,])

tL = 1/6
Length = 2
N = 30
fa = c()
fb = c()
fc = c()
ftime = c()
for(i in 1:(length(TimeGrid)-1)){
  timegrid = TimeGrid[c(i,(i+1))]
  TL1 = TL[c(i,(i+1))]
  TR1 = TR[c(i,(i+1))]
  PL1 = PL[c(i,(i+1))]
  PR1 = PR[c(i,(i+1))]
  res = wwwcost(N, tL, PL1, TL1, PR1, TR1, Length, timegrid)
  fa = c(fa,res[[1]][2:(length(res[[1]])-1)])
  fb = c(fb,res[[2]][2:(length(res[[2]])-1)])
  fc = c(fc,res[[3]][2:(length(res[[3]])-1)])
  ftime = c(ftime, res[[4]])

}
ad = as.matrix(data.frame(ftime=round(ftime,4),aa=round(fa,4), type=rep(1,length(ftime))))
bd = as.matrix(data.frame(ftime=round(ftime,4),bb=round(fb,4), type=rep(2,length(ftime))))
cd = as.matrix(data.frame(ftime=round(ftime,4),cc=round(fc,4), type=rep(3,length(ftime))))
plotdata = rbind(ad,bd)
plotdata = rbind(plotdata, cd)
plotdata = as.data.frame(plotdata)
plotdata[,3] = factor(plotdata[,3])
ggplot(plotdata, aes(x=ftime, y=aa, colour = type, shape = type)) +
  geom_point()# +facet_wrap( ~ type)
