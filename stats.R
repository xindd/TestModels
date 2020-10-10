

# 刺激信号的可能形式 h(t)

# 1、持续刺激信号，例：h(t) = k

# 2、逐渐衰减的信号，例：h(t) = exp(-t)

# 3、脉冲信号，例：h(t) =

# 4、其他形式的信号，例：h(t) = *

# 假设：在刺激之前速度处于稳定状态


# 拟合函数

# 刺激函数形式
#func = function(x){
#  exp(-x)
#}
#
## 根据阈值 m 筛选方差最小的拟合结果
#selectFit = function(x,y, num=100){
#  counters_m = spline(c(0,max(y)),c(1,2),num)$x
#  res = c()
#  for(m in counters_m){
#    tryCatch({
#      result=nls(y~f(x, a_s, S, m, func), start = list(a_s=0.1, S=0.1))
#      res = rbind(res, c(coefficients(result)['a_s'], coefficients(result)['S'], m, result$m$deviance()))
#    },error=function(e){break()}#{cat("ERROR :",conditionMessage(e),"\n")}
#    )
#  }
#  if( is.null(res)){
#    res=c()
#    minm=c(NA,NA,NA,NA)
#  }else{
#    colnames(res) = c('a_s','S','m','deviance')
#    minm = res[which(res[,'deviance']==min(res[,'deviance']))[1],]
#  }
#  
#  return(list(res,minm))
#}

# 刺激函数形式
signalFunc = function(t){
  h0=1.72290706
  h1= 0.66149998
  h2= 2.11254905
  t1= 1.10589967
  t2= 0.84831985
  bet=0.88761123
  #return(-bet*(h0 + (-h0 + h1)/(1 + exp(-bet*(t - t1))))*(h1 - h2)*exp(bet*(t - t2))/(h1*(exp(bet*(t - t2)) + 1)**2) + bet*(-h0 + h1)*(h2 + (h1 - h2)/(exp(bet*(t - t2)) + 1))*exp(-bet*(t - t1))/(h1*(1 + exp(-bet*(t - t1)))**2))
  #return(10*t/t)
  return(sin(t))
}


plot(x,signalFunc(x), type='l')
a.func<-function(t, parmaters, signalFunc){
  print(parmaters[2] * signalFunc(t) / (1 + exp(-parmaters[4]*(t-parmaters[3]))))
  parmaters[1] + parmaters[2] * signalFunc(t) / (1 + exp(-parmaters[4]*(t-parmaters[3])))
}


minE <- function(t,y,step=0.01,signalFunc, error=0.001,maxloop=100){
  time_len = length(t)
  # 初始化参数
  a_s = rnorm(1) 
  S = rnorm(1)
  t0 = rnorm(1)
  w = rnorm(1)
  xita = c(a_s, S, t0, w)
  param = xita
  energy = c()
  for(i in 1:maxloop){
    # 计算预测值
    ft = xita[1] + xita[2] * signalFunc(t) / (1 + exp(-xita[4]*(t-xita[3])))
    # 计算误差E
    E = 0.5*sum((y-ft)^2)
    energy = c(energy, E)
    energy_len = length(energy)
    
    # 计算偏导数
    du = c(
      sum(ft - y), # a_s
      sum((ft - y)*signalFunc(t)/(1 + exp(-xita[4]*(t-xita[3])))), # S
      sum((ft - y)*xita[2]*signalFunc(t)*xita[4]*exp(-xita[4]*(t-xita[3]))/(1 + exp(-xita[4]*(t-xita[3])))^2), # t0
      sum((ft - y)*xita[2]*signalFunc(t)*(xita[3]-t)*exp(-xita[4]*(t-xita[3]))/(1 + exp(-xita[4]*(t-xita[3])))^2) # w
    )
    # 校正参数
    if(energy_len>=2){
      if(abs(energy[energy_len]-energy[energy_len-1])<error #&
         #sum(abs(du)<error)==4
         ){
        print('break')
        break()
      }
    }
    #print(du)
    xita = xita - step*du
    param = rbind(param,xita)
  }
  return(list(param=param, energy=energy))
}

res = minE(x,y2,step=0.01,signalFunc, error=0,maxloop=2000)

par(mfrow=c(3,2))
plot(1:length(res$param[,1]), res$param[,1], ylab ='a_s')
plot(1:length(res$param[,2]), res$param[,2], ylab ='S')
plot(1:length(res$param[,3]), res$param[,3], ylab ='t0')
plot(1:length(res$param[,4]), res$param[,4], ylab ='w')
plot(res$energy)
print(res$param[length(res$param[,1]),])
plot(x,y2,type='l')
lines(x,a.func(x, res$param[length(res$param[,1]),], signalFunc))
plot(x,a.func(x, res$param[length(res$param[,1]),], signalFunc))

y = c(0.23113671, 0.10334868, 0.08452748, 0.29367640, 0.21612925)
y2 = c(4.705845,  3.970896,  6.928128, 12.676947, 12.112753)
y3 = c(4.722820, 3.923267, 4.189115, 6.171562, 5.811192)
y4 = c(12.407184, 15.481512, 18.224904, 10.284977,  9.447397)
x = c(1, 2, 4, 6, 8)

r1 = selectFit(x,y)

min_m = r1[[2]]

plot(x,y, type='l')
lines(x,f(x,min_m[1], min_m[2],min_m[3],func),type='l',col='red')


r2 = selectFit(x,y2)
min_m2 = r2[[2]]

plot(x,y2, type='l')
lines(x,f(x,min_m2[1], min_m2[2],min_m2[3],func),type='l',col='red')

#########################################################

# 过滤na行及方差过大的数据，波动异常值或无波动值
filter.outers <- function(data, fold=1.5){
  filter_na = na.omit(data)
  eval_var=apply(filter_na, 1, var)
  fenwei = quantile(eval_var)
  # 去除异常值：将个体是距离底线或顶线的距离超过fold倍的箱体高度视为离群值。
  filter_var <- filter_na[(eval_var < (fenwei[4] + (fenwei[4]-fenwei[2])*fold) & 
                             eval_var > (fenwei[2] - (fenwei[4]-fenwei[2])*fold)),]
  return(filter_var)
}

## filter zeros rates; filter rate < 0.5
filterzeros <- function(data, threshold=0.5){
  rownum = c()
  for(i in 1:dim(data)[1]){
    if(length(which(data[i,]==0)) >= 3){
      rownum = c(rownum, i)
    }else if(length(which(data[i,]<=threshold)) >=5){
      rownum = c(rownum, i)
    }
  }
  data=data[-rownum,]
}

at_filter_var = filter.outers(dlmat, fold=1.5)
at_filter_zero = filterzeros(at_filter_var)

boxplot(at_filter_zero)
summary(at_filter_zero)

at_result = t(apply(at_filter_zero,1,function(x){selectFit(t,x)[[2]]}))
pheatmap(na.omit(at_result[,2:3]),cluster_row=TRUE,
         cluster_col=FALSE,
         #scale="row",
         legend=TRUE,
         color = colorRampPalette(c("yellow", "red"))(100))



len = dim(dlmat)[1]
at_stats = c()
for(i in 1:len){
  r1 = selectFit(x,dlmat[i,])
  at_stats = rbind(at_stats, r1[[2]])
}
#t = c(0, 2, 4, 6, 8)
#da = data.frame(dlmat)
#result = apply(da,1,function(x){selectFit(t,x)[[2]]})
#summary(result)
#coefficients(result)
pheatmap(na.omit(at_result[,2:3]),cluster_row=TRUE,
         cluster_col=FALSE,
         scale="none",
         legend=TRUE,
         color = colorRampPalette(c("yellow", "red"))(100))

len = dim(dlmct)[1]
ct_stats = c()
for(i in 1:len){
  r1 = selectFit(x,dlmct[i,])
  ct_stats = rbind(ct_stats, r1[[2]])
}
#t = c(0, 2, 4, 6, 8)
#da = data.frame(dlmat)
result = apply(da[1:20,],1,function(x){selectFit(t,x)[[2]]})
#summary(result)
#coefficients(result)
pheatmap(na.omit(at_stats[,2:3]),cluster_row=TRUE,
         cluster_col=FALSE,
         scale="row",
         legend=TRUE,
         color = colorRampPalette(c("yellow", "red"))(100))

pheatmap(na.omit(ct_stats[,2:3]),cluster_row=FALSE,
         cluster_col=FALSE,
         scale="row",
         legend=TRUE,
         color = colorRampPalette(c("yellow", "red"))(100))

at_stats_filter = na.omit(at_stats[,2:3])

var_data = apply(at_stats_filter, 1, var)

sp <- boxplot(var_data)
at_stats_filter <- at_stats_filter[(at_stats_filter[,1] < (sp$stats[4,1] + (sp$stats[4,1]-sp$stats[2,1])*1) & 
                                              at_stats_filter[,1] > (sp$stats[2,1] - (sp$stats[4,1]-sp$stats[2,1])*1)),] # 去除异常值：将个体是距离底线或顶线的距离超过3倍的箱体高度视为离群值。

boxplot(at_stats_filter[,2])
