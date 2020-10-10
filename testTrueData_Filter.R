###############################################
library(dlm)
library(ggplot2)
library(grid)
library(reshape2)
################################

# error gene 9997, 9985, 9963, 9863
# gene 77595
TL = c(2.085375, 1.801374, 1.915773, 1.843418, 1.524426, 1.531284, 1.755409, 1.660484, 2.021102)
#TL = rpkms_rep1$foursu_exons['9985',]
TT = c(2.852587, 2.477313, 2.732977, 2.333752, 2.117182, 2.749777, 2.334962, 2.472333, 2.129490)
#TT = rpkms_rep1$total_exons['9985',]
PL = c(0.6361671, 0.6074711, 0.6403929, 0.6068481, 0.4185983, 0.4488358, 0.5062905, 0.5322454, 0.6000417)
#PL = rpkms_rep1$foursu_introns['9985',]
PT = c(0.4332703, 0.3704975, 0.3667127, 0.3228514, 0.2734701, 0.2819486, 0.3483869, 0.3103953, 0.2797586)
#PT = rpkms_rep1$total_introns['9985',]
at = c(7.528854, 4.980970, 5.876078, 5.658038, 4.855749, 4.002908, 4.661779, 3.858624, 5.182438)
ct = c(17.376805, 17.701023, 16.713811, 19.653820, 19.661656, 7.437069, 15.692097, 15.194234, 21.476798)
bt = c(3.111976, 4.7684460, 1.7874769, 4.085951, 3.156609, 2.0692682, 2.6770777, 2.2304268, 3.2218976)

################################
t_time = c(0, 0.17, 0.33, 0.5, 1, 2, 4, 8, 16)
tL <- 1/6

################################
trueData = data.frame(TL,TT,PL,PT)
trueData$id = t_time
colnames(trueData) = c('TL','TT','PL','PT','id')
trueData = melt(trueData, id.vars = 'id')

plot1 <- ggplot(trueData, aes(id, value)) + geom_line(aes(color=variable), size=1)

#inspectd = data.frame(at,ct,bt)
#inspectd$id = t_time
#inspectd = melt(inspectd, id.vars = 'id')

#plot2 <- ggplot(inspectd, aes(id, value)) + geom_line(aes(color=variable), size=1)
########
kalmanRates = function(t_time, tL, foursu_exons, total_exons, total_introns, splitnum=100, JW=0.0001){
  
  # 初始化模型
  #t_time2 = spline(t_time, 1:length(t_time), n = splitnum)$x
  num_data = length(t_time)
  dert = c(tL, diff(t_time))
  
  JGG <- matrix(c(0,1,0,0),nrow=2,byrow=TRUE)
  GG <- matrix(c(1,1,0,1),nrow=2,byrow=TRUE)
  FF <- matrix(c(1,0), nrow=1)
  m0 <- c(0,0)
  C0 <- diag(1,2)
  W <- diag(JW,2)
  V <- 1
  X = matrix(dert, nrow = num_data)
  X = cbind(X, rep(0,num_data))
  X = cbind(X, rep(1,num_data))
  my_dlm <- dlm(X=X, FF=FF,V=V, JGG=JGG, GG=GG, W=W,m0=m0, C0=C0)
  
  # 计算速率
  dlmrates = list()
  dlmres = list()
  foursu_exons = rbind(foursu_exons)
  total_exons = rbind(total_exons)
  total_introns = rbind(total_introns)
  
  for(gene in 1:dim(foursu_exons)[1]){
    TL = foursu_exons[gene,]#spline(t_time, foursu_exons[gene,], n = splitnum)$y
    TT = total_exons[gene,]#spline(t_time, total_exons[gene,], n = splitnum)$y
    PT = total_introns[gene,]#spline(t_time, total_introns[gene,], n = splitnum)$y
    
    pre_TT <- dlmSmooth(TT,my_dlm)
    pre_PT <- dlmSmooth(PT,my_dlm)
    
    dlmres=list(gene=list(pre_PT=dropFirst(pre_PT$s), pre_TT=dropFirst(pre_TT$s)))
    # 计算速率
    erP = dropFirst(pre_PT$s)[,1] / PT#PT*(PT-dropFirst(pre_PT$s)[,1])/(TT+PT)
    erT = dropFirst(pre_TT$s)[,1] / TT#TT*(TT-dropFirst(pre_TT$s)[,1])/(TT+PT)
    pre_TL = TL*(erP+erT)/2#TL-(erP+erT)/2
    at_dlm = pre_TL/tL
    # at_dlm = TL/tL
    ct_dlm = (at_dlm-dropFirst(pre_PT$s)[,2])/dropFirst(pre_PT$s)[,1]
    bt_dlm = (at_dlm-dropFirst(pre_TT$s)[,2])/(dropFirst(pre_TT$s)[,1]-dropFirst(pre_PT$s)[,1])
    dlmrates$at=rbind(dlmrates$at, at_dlm)
    dlmrates$ct=rbind(dlmrates$ct, ct_dlm)
    dlmrates$bt=rbind(dlmrates$bt, bt_dlm)
  }
  dlmrates$at[dlmrates$at<0] = NA
  dlmrates$ct[dlmrates$ct<0] = NA
  dlmrates$bt[dlmrates$bt<0] = NA
  rownames(dlmrates$at) = rownames(foursu_exons)
  rownames(dlmrates$ct) = rownames(foursu_exons)
  rownames(dlmrates$bt) = rownames(foursu_exons)
  return(list(rates=dlmrates,predata=dlmres,preTL=pre_TL))
  
}
kalmanRates2 = function(t_time, tL, foursu_exons, total_exons, total_introns, splitnum=100, JW=0.0001){
  
  # 初始化模型
  #t_time2 = spline(t_time, 1:length(t_time), n = splitnum)$x
  num_data = length(t_time)
  dert = c(tL, diff(t_time))
  
  JGG <- matrix(c(0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0),nrow=4,byrow=TRUE)
  GG <- matrix(c(1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,1),nrow=4,byrow=TRUE)
  
  FF <- matrix(c(1,0,0,0,0,0,1,0), nrow=2,byrow=TRUE)
  
  m0 <- c(0,0,0,0)
  C0 <- diag(1,4)
  
  W <- diag(JW,4)
  V <- diag(1,2)
  X = matrix(dert, nrow = num_data)
  X = cbind(X, rep(0,num_data))
  X = cbind(X, rep(1,num_data))
  my_dlm <- dlm(X=X, FF=FF,V=V, JGG=JGG, GG=GG, W=W,m0=m0, C0=C0)
  
  # 计算速率
  dlmrates = list()
  dlmres = list()
  foursu_exons = rbind(foursu_exons)
  total_exons = rbind(total_exons)
  total_introns = rbind(total_introns)
  for(gene in 1:dim(foursu_exons)[1]){
    TL = foursu_exons[gene,]#spline(t_time, foursu_exons[gene,], n = splitnum)$y
    TT = total_exons[gene,]#spline(t_time, total_exons[gene,], n = splitnum)$y
    PT = total_introns[gene,]#spline(t_time, total_introns[gene,], n = splitnum)$y
    
    pre <- dlmSmooth(cbind(PT,TT),my_dlm)
    
    dlmres=list(gene=dropFirst(pre$s))
    # 计算速率
    erP = dropFirst(pre$s)[,1] / PT#PT*(PT-dropFirst(pre_PT$s)[,1])/(TT+PT)
    erT = dropFirst(pre$s)[,3] / TT#TT*(TT-dropFirst(pre_TT$s)[,1])/(TT+PT)
    pre_TL = TL*(erP+erT)/2#TL-(erP+erT)/2
    at_dlm = pre_TL/tL
    # at_dlm = TL/tL
    ct_dlm = (at_dlm-dropFirst(pre$s)[,2])/dropFirst(pre$s)[,1]
    bt_dlm = (at_dlm-dropFirst(pre$s)[,4])/(dropFirst(pre$s)[,3]-dropFirst(pre$s)[,2])
    dlmrates$at=rbind(dlmrates$at, at_dlm)
    dlmrates$ct=rbind(dlmrates$ct, ct_dlm)
    dlmrates$bt=rbind(dlmrates$bt, bt_dlm)
  }
  dlmrates$at[dlmrates$at<0] = NA
  dlmrates$ct[dlmrates$ct<0] = NA
  dlmrates$bt[dlmrates$bt<0] = NA
  rownames(dlmrates$at) = rownames(foursu_exons)
  rownames(dlmrates$ct) = rownames(foursu_exons)
  rownames(dlmrates$bt) = rownames(foursu_exons)
  return(list(rates=dlmrates,predata=dlmres,preTL=pre_TL))
  
}

grid.newpage() # 生成新画布
pushViewport(viewport(layout = grid.layout(2,2))) # 指定一个3乘3的layout
print(plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) # 画
print(plot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) # 画
#################################

res = kalmanRates2(t_time, tL, TL, TT, PT, splitnum=100, JW=11)


dlmRes <- data.frame(res$predata$gene[,3], res$predata$gene[,4], res$predata$gene[1],res$predata$gene[2], res$preTL)
colnames(dlmRes) = c('pre_TT','pre_TT_dot', 'pre_PT', 'pre_PT_dot','pre_TL')
dlmRes$TT = TT #spline(t_time, TT, n = 100)$y
dlmRes$PT = PT #spline(t_time, PT, n = 100)$y
dlmRes$TL = TL #spline(t_time, TL, n = 100)$y
dlmRes$id <- t_time# spline(t_time, 1:length(t_time), n = 100)$x
dlmRes = melt(dlmRes, id.vars = 'id')
ggplot(dlmRes, aes(id, value)) + geom_line(aes(color=variable), size=1)


dlmRates = cbind(t(res$rates$at),t(res$rates$ct))
dlmRates = as.data.frame(cbind(dlmRates,t(res$rates$bt)))
colnames(dlmRates) = c('at','ct','bt')
dlmRates$id = t_time #spline(t_time, 1:length(t_time), n = 100)$x
dlmRates = melt(dlmRates, id.vars = 'id')
plot4 <- ggplot(dlmRates, aes(id, value)) + geom_line(aes(color=variable), size=1)
##################################################
print(plot3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) # 画
print(plot4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)) # 画
popViewport() # 退出当前位置











##################################################
# 方法二
true_a = TL/tL
truedata=matrix(c(PT,TT), nrow = num_data, byrow=FALSE)
Z = matrix(c(true_a,true_a), nrow = num_data, byrow=FALSE)
P_T_1 = c(0, PT[1:(num_data-1)])
T_T_1 = c(0, TT[1:(num_data-1)])
a_t_1 = c(0, true_a[1:(num_data-1)])

GG = diag(1,4)
JGG = diag(0,4)
JGG[1,1] = 1 
JGG[2,2] = 2

FF <- matrix(c(1,0,1,0,0,1,0,1),nrow=2,byrow=TRUE)
JFF = matrix(c(3,0,1,0,0,4,0,1),nrow=2,byrow=TRUE)

m0 <- c(0,0,0,0)
C0 <- diag(1,4)

w=.001
ww = cbind((true_a-a_t_1+w)/PT,(true_a-a_t_1+w)/(TT-PT),w,w)
WW = c()
for(i in 1:num_data){
  WW = rbind(WW,as.vector(as.matrix(ww[i,], nrow=4) %*% ww[i,]))
}
W <- diag(w,4)
JW = matrix(5:20, nrow = 4)
V <- diag(1,2)
X = matrix(c(P_T_1/PT, 
             (T_T_1-P_T_1)/(TT-PT),
             PT, 
             (TT-PT)),nrow=num_data,byrow=FALSE)
X = cbind(X, WW)

my_dlm <- dlm(X=X,JFF=FF,FF=FF,V=V, JGG=JGG, GG=GG, JW=JW,W=W,m0=m0, C0=C0)


y = dlmSmooth(Z,my_dlm)
pre_Z = matrix(c(0,0),nrow = 1)
for(k in 1:num_data){
  pre_Z=rbind(pre_Z,t(matrix(c(P_data[k],0,1,0,0,(T_data[k]-P_data[k]),0,1),nrow=2,byrow=TRUE)%*%as.matrix(y$s[k+1,],nrow=4)))
}
pre_Z = dropFirst(pre_Z)

df <- data.frame(truedata,pre_Z, Z, dropFirst(y$s))
colnames(df) = c('PT','TT', 'a1', 'a2', 'at1','at2', 'c','b','p_dot','t_dot')
df$id <- t_time
ggplot(data=df,aes(x=id,y=at1))+ # true_a
  geom_line()+
  geom_line(aes(y=at1),color="red")+ # a
  geom_line(aes(y=c),color="green")+ # c
  geom_line(aes(y=b),color="blue") # b





