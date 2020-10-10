

# 刺激信号的可能形式 h(t)

# 1、持续刺激信号，例：h(t) = k

# 2、逐渐衰减的信号，例：h(t) = exp(-t)

# 3、脉冲信号，例：h(t) =

# 4、其他形式的信号，例：h(t) = *

# 假设：在刺激之前速度处于稳定状态

#########################################
# 拟合函数
#f=function(t, a_s, S, m, func){
#  ifelse(func(t)>=m, a_s + S*func(t), a_s)
#}
#
## 刺激函数形式
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
#plot(x,signalFunc(x), type='l')
a.func<-function(t, parmaters, signalFunc){
  parmaters[1] + parmaters[2] * signalFunc(t) / (1 + exp(-(t-parmaters[3])))
}

# 估计参数 a_s，S，t
minE <- function(t,y,step=0.01,signalFun='down', error=0.001,maxloop=100, detail=FALSE){
  time_len = length(t)
  # 初始化参数
  a_s = rnorm(1) 
  S = rnorm(1)
  t0 = rnorm(1)
  xita = c(a_s, S, t0)
  if(detail){
    param = xita
    energy = c()
  }
  if(signalFun=='down'){
    # 信号函数形式
    signalFunc = function(t){
      h0=1.72290706
      h1= 0.66149998
      h2= 2.11254905
      t1= 1.10589967
      t2= 0.84831985
      bet=0.88761123
      #return(-bet*(h0 + (-h0 + h1)/(1 + exp(-bet*(t - t1))))*(h1 - h2)*exp(bet*(t - t2))/(h1*(exp(bet*(t - t2)) + 1)**2) + bet*(-h0 + h1)*(h2 + (h1 - h2)/(exp(bet*(t - t2)) + 1))*exp(-bet*(t - t1))/(h1*(1 + exp(-bet*(t - t1)))**2))
      #return(10*t/t)
      #return(exp(-t))
      return(exp(-x))
    }
  }
  
  for(i in 1:maxloop){
    # 计算预测值
    ft = xita[1] + xita[2] * signalFunc(t) / (1 + exp(-(t-xita[3])))
    # 计算误差E
    E = 0.5*sum((y-ft)^2)
    if(detail){
      energy = c(energy, E)
    }
    
    #energy_len = length(energy)
    
    # 计算偏导数
    du = c(
      sum(ft - y), # a_s
      sum((ft - y)*signalFunc(t)/(1 + exp(-(t-xita[3])))), # S
      sum((ft - y)*xita[2]*signalFunc(t)*exp(-(t-xita[3]))/(1 + exp(-(t-xita[3])))^2) # t0
    )
    # 校正参数
    #print(energy)
    #if(energy_len>=2){
    #  if(abs(energy[energy_len]-energy[energy_len-1])<error #&
    #     #sum(abs(du)<error)==4
    #  ){
    #    print('break')
    #    break()
    #  }
    #}
    #print(du)
    xita = xita - step*du
    if(detail){
      param = rbind(param,xita)
    }
   
  }
  if(detail){
    return(list(param=param, energy=energy))
  }else{
    return(xita)
  }
  
}


##########################################################
#测试
#res = minE(x,y3,step=0.1,signalFunc, error=0,maxloop=1000)
#par(mfrow=c(3,2))
#plot(1:length(res$param[,1]), res$param[,1], ylab ='a_s')
#plot(1:length(res$param[,2]), res$param[,2], ylab ='S')
#plot(1:length(res$param[,3]), res$param[,3], ylab ='t0')
#plot(res$energy)
#print(res$param[length(res$param[,1]),])
#plot(x,y3,type='l')
#lines(x,a.func(x, res$param[length(res$param[,1]),], signalFunc))
#plot(x,a.func(x, res$param[length(res$param[,1]),], signalFunc))
#y = c(0.23113671, 0.10334868, 0.08452748, 0.29367640, 0.21612925)
#y2 = c(4.705845,  3.970896,  6.928128, 12.676947, 12.112753)
#y3 = c(4.722820, 4.923267, 5.189115, 6.171562, 5.811192)
#y4 = c(18.407184, 15.481512, 13.224904, 10.284977,  9.447397)
#x = c(1, 2, 4, 6, 8)
#
#r1 = selectFit(x,y)
#
#min_m = r1[[2]]
#
#plot(x,y, type='l')
#lines(x,f(x,min_m[1], min_m[2],min_m[3],func),type='l',col='red')
#
#
#r2 = selectFit(x,y2)
#min_m2 = r2[[2]]
#
#plot(x,y2, type='l')
#lines(x,f(x,min_m2[1], min_m2[2],min_m2[3],func),type='l',col='red')

#########################################################

# 过滤速率数据
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
    if(length(which(data[i,]==0)) >= ceiling(dim(data)[2]/2)){
      rownum = c(rownum, i) # 半数时间点均为0表达值去掉
    }else if(length(which(data[i,]<=threshold)) >=dim(data)[2]){
      rownum = c(rownum, i) # 表达值均小于threshold的去掉
    }
  }
  data=data[-rownum,]
}


testfun=function(y,t, step=0.1,signalFun='down', error=0,maxloop=1000){
  return(minE(t,y,step=step,signalFun=signalFun, error=error,maxloop=maxloop))
  
}

######################
# 处理速率数据
handle.rates <- function(ttime,dlmat, dlmct, dlmbt){
  # 转录速率at
  message('正在计算转录at速率...')
  at_filter_var = filter.outers(dlmat, fold=1.5)
  at_filter_zero = filterzeros(at_filter_var)
  
  at_as = t(apply(at_filter_zero, 1, testfun, t=ttime, step=0.1))
  
  #par(mfrow=c(1,1))
  #boxplot(at_filter_zero)
  #boxplot(at_as)
  #
  #pheatmap(na.omit(at_as[,1:3]),cluster_row=FALSE,
  #         cluster_col=FALSE,
  #         scale="row",
  #         legend=TRUE,
  #         color = colorRampPalette(c("yellow", "red"))(100))
  
  #########################################
  # 加工速率ct
  message('正在计算加工ct速率...')
  ct_filter_var = filter.outers(dlmct, fold=1.5)
  ct_filter_zero = filterzeros(ct_filter_var)
  #boxplot(ct_filter_zero)
  
  ct_cs = t(apply(ct_filter_zero, 1, testfun, t=ttime, step=0.1))
  
  #boxplot(ct_cs)
  #
  #pheatmap(na.omit(ct_cs[,1:3]),cluster_row=TRUE,
  #         cluster_col=FALSE,
  #         scale="row",
  #         legend=TRUE,
  #         color = colorRampPalette(c("yellow", "red"))(100))
  #
  ##########################################
  # 降解速率bt
  message('正在计算降解bt速率...')
  bt_filter_var = filter.outers(dlmbt, fold=1.5)
  bt_filter_zero = filterzeros(bt_filter_var)
  #boxplot(bt_filter_zero)
  
  bt_bs = t(apply(bt_filter_zero, 1, testfun, t=t_time, step=0.1))
  #boxplot(bt_bs)
  #
  #pheatmap(na.omit(bt_bs[,1:3]),cluster_row=TRUE,
  #         cluster_col=FALSE,
  #         scale="row",
  #         legend=TRUE,
  #         color = colorRampPalette(c("yellow", "red"))(100))
  
  ###########################################
  
  # 基因取交集
  ratesgenelist = intersect(rownames(at_as), rownames(ct_cs))
  ratesgenelist = intersect(ratesgenelist, rownames(bt_bs))
  # 合并稳态速率
  message('合并稳态速率...')
  stats_rates = cbind(at_as[ratesgenelist,1], ct_cs[ratesgenelist,1])
  stats_rates = cbind(stats_rates, bt_bs[ratesgenelist,1])
  # 合并敏感性
  sensitive = cbind(at_as[ratesgenelist,2], ct_cs[ratesgenelist,2])
  sensitive = cbind(sensitive, bt_bs[ratesgenelist,2])
  
  #boxplot(stats_rates)
  #
  #pheatmap(na.omit(stats_rates),cluster_row=TRUE,
  #         cluster_col=FALSE,
  #         scale="row",
  #         legend=TRUE,
  #         color = colorRampPalette(c("yellow", "red"))(100))
  #
  #########################################
  # 稳态基因表达值
  message('计算稳态基因表达值...')
  stats_expression = function(x){
    if(x[2]!=0 & x[3]!=0){
      total_p_s = x[1]/x[2]
      total_t_s = x[1]/x[3] + total_p_s
      return(c(total_p_s, total_t_s))
    }else{
      return(c(0,0))
    }
  }
  
  base_expression = t(apply(stats_rates,1,stats_expression))
  
  base_gene = rownames(base_expression)
  
  observer_p = totintr_inspect[base_gene,]
  observer_t = totexon_inspect[base_gene,]
  
  # 计算（基因观测表达值-基础表达值）的方差， 1.5倍四分位外的为差异表达
  
  message('筛选差异表达基因...')
  gene_var = apply(observer_p-base_expression[,1],1,var)
  qu = quantile(gene_var)
  up_gene_p = observer_p[names(which((gene_var > (qu[4] + (qu[4]-qu[2])*1.5))==TRUE)),]
  down_gene_P = observer_p[names(which((gene_var < (qu[2] - (qu[4]-qu[2])*1.5))==TRUE)),]
  
  gene_var = apply(observer_t-base_expression[,2],1,var)
  # plot(gene_var)
  qu = quantile(gene_var)
  up_gene_t = observer_t[names(which((gene_var > (qu[4] + (qu[4]-qu[2])*1.5))==TRUE)),]
  down_gene_t = observer_t[names(which((gene_var < (qu[2] - (qu[4]-qu[2])*1.5))==TRUE)),]
  
  # 提取在 p，t均差异表达的基因
  filter_gene_final = intersect(rownames(up_gene_p), rownames(up_gene_t))
  filter_gene_final_t = observer_t[filter_gene_final,]
  filter_gene_final_p = observer_p[filter_gene_final,]
  
  return(list(t=filter_gene_final_t, p=filter_gene_final_p, rate=stats_rates,sensitive=sensitive))
  #pheatmap(na.omit(filter_gene_final_p),cluster_row=TRUE,
  #         cluster_col=FALSE,
  #         scale="row",
  #         legend=TRUE,
  #         color = colorRampPalette(c("yellow", "red"))(100))
  
}


# 分析rep相关性
ttime=c(0,2,4,6,8)
res1=handle.rates(ttime,dlmat, dlmct, dlmbt)
res2=handle.rates(ttime,dlmat2, dlmct2, dlmbt2)

# 基因overlap
overlap_gene = intersect(rownames(res1$t),rownames(res2$t))
library(grid)
library(futile.logger)
library(VennDiagram)
venn.diagram(x=list(rep1=rownames(res1$t), rep2=rownames(res2$t)),filename = "My3.png")
# 速率相关性
par(mfrow=c(2,3))
plot(res1$rate[overlap_gene,1], res2$rate[overlap_gene,1], ylab ='a2', xlab='a1')
plot(res1$rate[overlap_gene,2], res2$rate[overlap_gene,2], ylab ='c2', xlab='c1')
plot(res1$rate[overlap_gene,3], res2$rate[overlap_gene,3], ylab ='b2', xlab='b1')

plot(res1$sensitive[overlap_gene,1], res2$sensitive[overlap_gene,1], ylab ='Sa2', xlab='Sa1')
plot(res1$sensitive[overlap_gene,2], res2$sensitive[overlap_gene,2], ylab ='Sc2', xlab='Sc1')
plot(res1$sensitive[overlap_gene,3], res2$sensitive[overlap_gene,3], ylab ='Sb2', xlab='Sb1')


cor(res1$sensitive[overlap_gene,1], res2$sensitive[overlap_gene,1])
cor(res1$sensitive[overlap_gene,2], res2$sensitive[overlap_gene,2])
cor(res1$sensitive[overlap_gene,3], res2$sensitive[overlap_gene,3])


par(mfrow=c(1,1))

at_filter_var = filter.outers(dlmat, fold=1.5)
at_filter_zero = filterzeros(at_filter_var)

ct_filter_var = filter.outers(dlmct, fold=1.5)
ct_filter_zero = filterzeros(ct_filter_var)

bt_filter_var = filter.outers(dlmbt, fold=1.5)
bt_filter_zero = filterzeros(bt_filter_var)
boxplot(bt_filter_zero)
