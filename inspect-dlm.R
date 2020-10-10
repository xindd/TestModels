###############################################
library(dlm)
library(ggplot2)
library(grid)
library(reshape2)
################################
library(BiocParallel)
library(Biobase)
library(BiocGenerics)
library(parallel)
library(INSPEcT)
############读取表达数据####################
# 读取表达文件

rpkms_rep1 = list(foursu_exons=c(),
                  total_exons=c(),
                  foursu_introns=c(),
                  total_introns=c()
)
rpkms_rep2 = list(foursu_exons=c(),
                  total_exons=c(),
                  foursu_introns=c(),
                  total_introns=c()
)
dirpath = 'D:/output'
filename = list.files(dirpath)
rep1_geneid=c()
rep2_geneid=c()
for(name in filename){
  if(substr(name,12,12)=="1"){
    tmpfile = read.table(paste(dirpath,name,sep='/'), header = TRUE)
    tmpfile = tmpfile[order(tmpfile$id),]
    rpkms_rep1$foursu_exons=cbind(rpkms_rep1$foursu_exons,tmpfile$foursu_exons)
    rpkms_rep1$foursu_introns=cbind(rpkms_rep1$foursu_introns,tmpfile$foursu_introns)
    rpkms_rep1$total_exons=cbind(rpkms_rep1$total_exons,tmpfile$total_exons)
    rpkms_rep1$total_introns=cbind(rpkms_rep1$total_introns,tmpfile$total_introns)
    
  }
  if(substr(name,12,12)=="2"){
    tmpfile = read.table(paste(dirpath,name,sep='/'), header = TRUE)
    tmpfile = tmpfile[order(tmpfile$id),]
    rpkms_rep2$foursu_exons=cbind(rpkms_rep2$foursu_exons,tmpfile$foursu_exons)
    rpkms_rep2$foursu_introns=cbind(rpkms_rep2$foursu_introns,tmpfile$foursu_introns)
    rpkms_rep2$total_exons=cbind(rpkms_rep2$total_exons,tmpfile$total_exons)
    rpkms_rep2$total_introns=cbind(rpkms_rep2$total_introns,tmpfile$total_introns)
  }
}
rownames(rpkms_rep1$foursu_exons) = tmpfile$id
rownames(rpkms_rep1$foursu_introns) = tmpfile$id
rownames(rpkms_rep1$total_exons) = tmpfile$id
rownames(rpkms_rep1$total_introns) = tmpfile$id
rownames(rpkms_rep2$foursu_exons) = tmpfile$id
rownames(rpkms_rep2$foursu_introns) = tmpfile$id
rownames(rpkms_rep2$total_exons) = tmpfile$id
rownames(rpkms_rep2$total_introns) = tmpfile$id

colnames(rpkms_rep1$foursu_exons) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep1$foursu_introns) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep1$total_exons) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')
colnames(rpkms_rep1$total_introns) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')
colnames(rpkms_rep2$foursu_exons) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep2$foursu_introns) = c('4su_t_0','4su_t_2','4su_t_4','4su_t_6','4su_t_8')
colnames(rpkms_rep2$total_exons) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')
colnames(rpkms_rep2$total_introns) = c('total_t_0','total_t_2','total_t_4','total_t_6','total_t_8')

############提取部分基因####################
## filter zeros expression; filter expression < 0.5
filterzeros <- function(data){
  rownum = c()
  for(i in 1:dim(data)[1]){
    if(length(which(data[i,]==0)) >= 5){
      rownum = c(rownum, i)
    }else if(length(which(data[i,]<=0.1)) >=5){
      rownum = c(rownum, i)
    }
  }
  data=data[-rownum,]
}

# fielter error data , must have TL<TT, PT<TT
filtererror <- function(TL, TT, PT){
  w = dim(TL)[2]
  h = dim(TL)[1]
  rownum = c()
  for(i in i:h){
    if(sum(TL[i,]<TT[i,]) < w){
      rownum = c(rownum, i)
    }else if(sum(PT[i,]<TT[i,]) < w){
      rownum = c(rownum, i)
    }
  }
  data = list()
  data$tl = TL[-rownum, ]
  data$tt = TT[-rownum, ]
  data$pt = PT[-rownum, ]
  return(data)
}

labexon = rpkms_rep1$foursu_exons
labintr = rpkms_rep1$foursu_introns
totexon = rpkms_rep1$total_exons
totintr = rpkms_rep1$total_introns

labexon = filterzeros(labexon)
labintr = filterzeros(labintr)
totexon = filterzeros(totexon)
totintr = filterzeros(totintr)

dim(labexon)

#labexon = labexon[order(labexon[,1],labexon[,2],
#                        labexon[,3],labexon[,4],labexon[,5],decreasing = TRUE),][10:3000,]
#labintr = labintr[order(labintr[,1],labintr[,2],
#                        labintr[,3],labintr[,4],labintr[,5],decreasing = TRUE),][10:3000,]
#totexon = totexon[order(totexon[,1],totexon[,2],
#                        totexon[,3],totexon[,4],totexon[,5],decreasing = TRUE),][10:3000,]
#totintr = totintr[order(totintr[,1],totintr[,2],
#                        totintr[,3],totintr[,4],totintr[,5],decreasing = TRUE),][10:3000,]


labexon2 = rpkms_rep2$foursu_exons
labintr2 = rpkms_rep2$foursu_introns
totexon2 = rpkms_rep2$total_exons
totintr2 = rpkms_rep2$total_introns

labexon2 = filterzeros(labexon2)
labintr2 = filterzeros(labintr2)
totexon2 = filterzeros(totexon2)
totintr2 = filterzeros(totintr2)

#labexon2 = labexon2[order(labexon2[,1],labexon2[,2],
#                          labexon2[,3],labexon2[,4],labexon2[,5],decreasing = TRUE),][10:3000,]
#labintr2 = labintr2[order(labintr2[,1],labintr2[,2],
#                          labintr2[,3],labintr2[,4],labintr2[,5],decreasing = TRUE),][10:3000,]
#totexon2 = totexon2[order(totexon2[,1],totexon2[,2],
#                          totexon2[,3],totexon2[,4],totexon2[,5],decreasing = TRUE),][10:3000,]
#totintr2 = totintr2[order(totintr2[,1],totintr2[,2],
#                          totintr2[,3],totintr2[,4],totintr2[,5],decreasing = TRUE),][10:3000,]


genelist2 = intersect(rownames(labexon2),rownames(labintr2))
genelist2 = intersect(genelist2,rownames(totexon2))
genelist2 = intersect(genelist2,rownames(totintr2))


genelist = intersect(rownames(labexon),rownames(labintr))
genelist = intersect(genelist,rownames(totexon))
genelist = intersect(genelist,rownames(totintr))

genelist = intersect(genelist, genelist2)

labexon = labexon[genelist,]
labintr = labintr[genelist,]
totexon = totexon[genelist,]
totintr = totintr[genelist,]
labexon = labexon[order(rownames(labexon),decreasing = TRUE),]
labintr = labintr[order(rownames(labintr),decreasing = TRUE),]
totexon = totexon[order(rownames(totexon),decreasing = TRUE),]
totintr = totintr[order(rownames(totintr),decreasing = TRUE),]

labexon2 = labexon2[genelist,]
labintr2 = labintr2[genelist,]
totexon2 = totexon2[genelist,]
totintr2 = totintr2[genelist,]
labexon2 = labexon2[order(rownames(labexon2),decreasing = TRUE),]
labintr2 = labintr2[order(rownames(labintr2),decreasing = TRUE),]
totexon2 = totexon2[order(rownames(totexon2),decreasing = TRUE),]
totintr2 = totintr2[order(rownames(totintr2),decreasing = TRUE),]


#res1 = filtererror(labexon, totexon, totintr)
#labexon = res1$tl
#totexon = res1$tt
#totintr = res1$pt
#
#res2 = filtererror(labexon2, totexon2, totintr2)
#labexon2 = res2$tl
#totexon2 = res2$tt
#totintr2 = res2$pt
#
#
#genelist2 = intersect(rownames(labexon2),rownames(labintr2))
#genelist2 = intersect(genelist2,rownames(totexon2))
#genelist2 = intersect(genelist2,rownames(totintr2))


#genelist = intersect(rownames(labexon),rownames(labintr))
#genelist = intersect(genelist,rownames(totexon))
#genelist = intersect(genelist,rownames(totintr))
#
#genelist = intersect(genelist, genelist2)
#
#labexon = labexon[genelist,]
#labintr = labintr[genelist,]
#totexon = totexon[genelist,]
#totintr = totintr[genelist,]
#labexon = labexon[order(rownames(labexon),decreasing = TRUE),]
#labintr = labintr[order(rownames(labintr),decreasing = TRUE),]
#totexon = totexon[order(rownames(totexon),decreasing = TRUE),]
#totintr = totintr[order(rownames(totintr),decreasing = TRUE),]
#
#labexon2 = labexon2[genelist,]
#labintr2 = labintr2[genelist,]
#totexon2 = totexon2[genelist,]
#totintr2 = totintr2[genelist,]
#labexon2 = labexon2[order(rownames(labexon2),decreasing = TRUE),]
#labintr2 = labintr2[order(rownames(labintr2),decreasing = TRUE),]
#totexon2 = totexon2[order(rownames(totexon2),decreasing = TRUE),]
#totintr2 = totintr2[order(rownames(totintr2),decreasing = TRUE),]

############INSPECT速率###################
t_time <- c(0, 2, 4, 6, 8)
tL <- 1
num_data = length(t_time)
dert = c(0, diff(t_time))

myrates <- newINSPEcT(t_time, tL, labexon+labintr, totexon+totintr,
                      labintr, totintr, BPPARAM=SerialParam())
myrates2 <- newINSPEcT(t_time, tL, labexon2+labintr2, totexon2+totintr2,
                       labintr2, totintr2, BPPARAM=SerialParam())

##############dlm 参数初始化####################
#JGG <- matrix(c(0,1,0,0),nrow=2,byrow=TRUE)
#GG <- matrix(c(1,1,0,1),nrow=2,byrow=TRUE)
#FF <- matrix(c(1,0), nrow=1)
#m0 <- c(0,0)
#C0 <- diag(1,2)
#W <- diag(0.0001,2)
#V <- 1
#X = matrix(dert, nrow = num_data)
#X = cbind(X, rep(0,num_data))
#X = cbind(X, rep(1,num_data))
#my_dlm <- dlm(X=X, FF=FF,V=V, JGG=JGG, GG=GG, W=W,m0=m0, C0=C0)
##############dlm 速率计算#####
#dlmRates = function(t_time, tL, model, foursu_exons, total_exons, total_introns){
#  num_data = length(t_time)
#  dert = c(0, diff(t_time))
#  dlmrates = list()
#  for(gene in 1:dim(foursu_exons)[1]){
#    TL = foursu_exons[gene,]
#    TT = total_exons[gene,]
#    PT = total_introns[gene,]
#    
#    pre_TT <- dlmSmooth(TT,model)
#    pre_PT <- dlmSmooth(PT,model)
#    # 计算速率
#    at_dlm = TL/tL
#    ct_dlm = (at_dlm-dropFirst(pre_PT$s)[,2])/dropFirst(pre_PT$s)[,1]
#    bt_dlm = (at_dlm-dropFirst(pre_TT$s)[,2])/(dropFirst(pre_TT$s)[,1]-dropFirst(pre_PT$s)[,1])
#    dlmrates$at=rbind(dlmrates$at, at_dlm)
#    dlmrates$ct=rbind(dlmrates$ct, ct_dlm)
#    dlmrates$bt=rbind(dlmrates$bt, bt_dlm)
#  }
#  dlmrates$at[dlmrates$at<0] = NA
#  dlmrates$ct[dlmrates$ct<0] = NA
#  dlmrates$bt[dlmrates$bt<0] = NA
#  rownames(dlmrates$at) = rownames(foursu_exons)
#  rownames(dlmrates$ct) = rownames(foursu_exons)
#  rownames(dlmrates$bt) = rownames(foursu_exons)
#  return(dlmrates)
#}
t_time <- c(0, 2, 4, 6, 8)
tL <- 1
ins_gene = rownames(ratesFirstGuess(myrates, 'synthesis'))
labexon_inspect = labexon[ins_gene,]
totexon_inspect = totexon[ins_gene,]
totintr_inspect = totintr[ins_gene,]

labexon_inspect2 = labexon2[ins_gene,]
totexon_inspect2 = totexon2[ins_gene,]
totintr_inspect2 = totintr2[ins_gene,]

#dlmrates <- dlmRates(t_time, tL, my_dlm, labexon_inspect, totexon_inspect, totintr_inspect)
#dlmrates2 <- dlmRates(t_time, tL, my_dlm, labexon_inspect2, totexon_inspect2, totintr_inspect2)

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
  return(dlmrates)#list(rates=dlmrates,predata=dlmres,preTL=pre_TL))
  
}

kalmanRates2 = function(t_time, tL, foursu_exons, total_exons, total_introns, splitnum=100, JW=0.0001, JV=1){
  
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
  V <- diag(JV,2)
  X = matrix(dert, nrow = num_data)
  X = cbind(X, rep(0,num_data))
  X = cbind(X, rep(1,num_data))
  my_dlm <- dlm(X=X, FF=FF,V=V, JGG=JGG, GG=GG, W=W,m0=m0, C0=C0)
  
  # 计算速率
  dlmrates = list()
  #dlmres = list()
  foursu_exons = rbind(foursu_exons)
  total_exons = rbind(total_exons)
  total_introns = rbind(total_introns)
  for(gene in 1:dim(foursu_exons)[1]){
    TL = foursu_exons[gene,]#spline(t_time, foursu_exons[gene,], n = splitnum)$y
    TT = total_exons[gene,]#spline(t_time, total_exons[gene,], n = splitnum)$y
    PT = total_introns[gene,]#spline(t_time, total_introns[gene,], n = splitnum)$y
    
    pre <- dlmSmooth(cbind(PT,TT),my_dlm)
    
    #dlmres=list(gene=dropFirst(pre$s))
    # 计算速率
    erP = dropFirst(pre$s)[,1] / PT#PT*(PT-dropFirst(pre_PT$s)[,1])/(TT+PT)
    erT = dropFirst(pre$s)[,3] / TT#TT*(TT-dropFirst(pre_TT$s)[,1])/(TT+PT)
    pre_TL = TL*(erP+erT)/2#TL-(erP+erT)/2
    at_dlm = pre_TL/tL
    # at_dlm = TL/tL
    ct_dlm = (at_dlm-dropFirst(pre$s)[,2])/dropFirst(pre$s)[,1]
    bt_dlm = (at_dlm-dropFirst(pre$s)[,4])/(dropFirst(pre$s)[,3]-dropFirst(pre$s)[,1])
    dlmrates$at=rbind(dlmrates$at, at_dlm)
    dlmrates$ct=rbind(dlmrates$ct, ct_dlm)
    dlmrates$bt=rbind(dlmrates$bt, bt_dlm)
  }
  dlmrates$at[dlmrates$at<0] = NA
  dlmrates$ct[dlmrates$ct<0] = NA
  dlmrates$bt[dlmrates$bt<0] = NA
  
  dlmrates$at[dlmrates$at==Inf]=NA
  dlmrates$ct[dlmrates$ct==Inf]=NA
  dlmrates$bt[dlmrates$bt==Inf]=NA
  rownames(dlmrates$at) = rownames(foursu_exons)
  rownames(dlmrates$ct) = rownames(foursu_exons)
  rownames(dlmrates$bt) = rownames(foursu_exons)
  return(dlmrates)
  
}


##########选择最优方差#############
####### 提取高表达10%基因###### ####
sp <- boxplot(labexon_inspect)
labexon_inspect <- labexon_inspect[rownames(labexon_inspect < (sp$stats[4] + (sp$stats[4]-sp$stats[2])*3) & 
                                     labexon_inspect > (sp$stats[2] - (sp$stats[4]-sp$stats[2])*3)),] # 去除异常值：将个体是距离底线或顶线的距离超过3倍的箱体高度视为离群值。


labexon_inspect = labexon_inspect[order(labexon_inspect[,1],labexon_inspect[,2],
                                        labexon_inspect[,3],labexon_inspect[,4],labexon_inspect[,5],decreasing = TRUE),][100:200,]

genename = rownames(labexon_inspect)
totexon_inspect = totexon_inspect[genename,]
totintr_inspect = totintr_inspect[genename,]

labexon_inspect2 = labexon_inspect2[genename,]
totexon_inspect2 = totexon_inspect2[genename,]
totintr_inspect2 = totintr_inspect2[genename,]

corrates=c()
col_k=5
for(w in c(0,0.0001, 0.001, 0.01, 0.1, 0.5, 0.8, 1,1:20)){
  res1 = kalmanRates2(t_time, tL, labexon_inspect, totexon_inspect, totintr_inspect, splitnum=100, JW=w , JV=1)
  res2 = kalmanRates2(t_time, tL, labexon_inspect2, totexon_inspect2, totintr_inspect2, splitnum=100, JW=w , JV=1)
  # dlm
  dlmat = res1$at
  dlmat2 = res2$at
  dlmct = res1$ct
  dlmct2 = res2$ct
  dlmbt = res1$bt
  dlmbt2 = res2$bt
  
  ls = na.omit(data.frame(b1=dlmbt[,col_k], b2=dlmbt2[,col_k]))
  bcor = cor(ls[,1], ls[,2])
  
  ls = na.omit(data.frame(b1=dlmct[,col_k], b2=dlmct2[,col_k]))
  ccor = cor(ls[,1], ls[,2])
  
  ls = na.omit(data.frame(b1=dlmat[,col_k], b2=dlmat2[,col_k]))
  acor = cor(ls[,1], ls[,2])
  
  
  corrates = cbind(corrates, c(acor, bcor, ccor))
}

plot(corrates[1,], type='l', ylim = c(0,1), col = "red")
lines(corrates[2,],col = "green")
lines(corrates[3,],col = "black")
#library(scatterplot3d)
#my.mat=corrates
#s3d.dat=data.frame(columns=c(row(my.mat)),rows=c(col(my.mat)),value=c(my.mat))
#scatterplot3d(s3d.dat,type="h",lwd=5,pch="",x.ticklabs=colnames(my.mat),
#              y.ticklabs=rownames(my.mat), main="3D barplot")



res1 = kalmanRates2(t_time, tL, labexon_inspect, totexon_inspect+totintr_inspect, totintr_inspect, splitnum=100, JW=0.01)
res2 = kalmanRates2(t_time, tL, labexon_inspect2, totexon_inspect2+totintr_inspect2, totintr_inspect2, splitnum=100, JW=0.01)

# dlm
dlmat = res1$at
dlmat2 = res2$at
#plot(dlmat[,1],dlmat2[,1])
dlmct = res1$ct
dlmct2 = res2$ct
#plot(dlmct[,2],dlmct2[,2])
dlmbt = res1$bt
dlmbt2 = res2$bt
#plot(dlmbt[,1],dlmbt2[,1])
#cor(dlmbt[,1],dlmbt2[,1])

#横向查看基因相关性
#tmp = function(x){
#  cor(x[1:5],x[6:10])
#}
#at_cor = apply(cbind(dlmat,dlmat2), 1,tmp)
#ct_cor = apply(cbind(dlmct,dlmct2), 1,tmp)
#bt_cor = apply(cbind(dlmbt,dlmbt2), 1,tmp)
#plot(at_cor)
#plot(ct_cor)
#plot(bt_cor)

# inspect
synthe = ratesFirstGuess(myrates, 'synthesis')
synthe2 = ratesFirstGuess(myrates2, 'synthesis')

progress = ratesFirstGuess(myrates, 'processing')
progress2 = ratesFirstGuess(myrates2, 'processing')

degradation = ratesFirstGuess(myrates, 'degradation')
degradation2 = ratesFirstGuess(myrates2, 'degradation')

#synthe_cor = apply(cbind(synthe,synthe2), 1,tmp)
#progress_cor = apply(cbind(progress,progress2), 1,tmp)
#degradation_cor = apply(cbind(degradation,degradation2), 1,tmp)
#plot(synthe_cor)
#plot(progress_cor)
#plot(degradation_cor)

gl = intersect(rownames(synthe),rownames(synthe2))
gl = intersect(gl,rownames(progress))
gl = intersect(gl,rownames(progress2))
gl = intersect(gl,rownames(degradation))
gl = intersect(gl,rownames(degradation2))
synthe = synthe[gl,]
synthe2 = synthe2[gl,]
progress = progress[gl,]
progress2 = progress2[gl,]
degradation = degradation[gl,]
degradation2 = degradation2[gl,]
#排序
synthe = synthe[order(rownames(synthe),decreasing = TRUE),]
synthe2 = synthe2[order(rownames(synthe2),decreasing = TRUE),]
progress = progress[order(rownames(progress),decreasing = TRUE),]
progress2 = progress2[order(rownames(progress2),decreasing = TRUE),]
degradation = degradation[order(rownames(degradation),decreasing = TRUE),]
degradation2 = degradation2[order(rownames(degradation2),decreasing = TRUE),]

#plot(synthe[,1],synthe2[,1])
#plot(progress[,1],progress2[,1])
#plot(degradation[,1],degradation2[,1])
################################

col_k = 2

plot1 <- ggplot(na.omit(data.frame(a1=dlmat[,col_k], a2=dlmat2[,col_k])), aes(a1, a2)) + geom_point()
plot2 <- ggplot(na.omit(data.frame(c1=dlmct[,col_k], c2=dlmct2[,col_k])), aes(c1, c2)) + geom_point()# +xlim(0,2000)+ylim(0,2000)
plot3 <- ggplot(na.omit(data.frame(b1=dlmbt[,col_k], b2=dlmbt2[,col_k])), aes(b1, b2)) + geom_point()# +xlim(0,125)+ylim(0,100)

plot5 <- ggplot(na.omit(data.frame(a1=synthe[,col_k], a2=synthe2[,col_k])), aes(a1, a2)) + geom_point()
plot6 <- ggplot(na.omit(data.frame(c1=progress[,col_k], c2=progress2[,col_k])), aes(c1, c2)) + geom_point()
plot7 <- ggplot(na.omit(data.frame(b1=degradation[,col_k], b2=degradation2[,col_k])), aes(b1, b2)) + geom_point()


ls = na.omit(data.frame(b1=dlmbt[,col_k], b2=dlmbt2[,col_k]))
cor(ls[,1], ls[,2])
ls = na.omit(data.frame(b1=degradation[,col_k], b2=degradation2[,col_k]))
cor(ls[,1], ls[,2])

ls = na.omit(data.frame(b1=dlmct[,col_k], b2=dlmct2[,col_k]))
cor(ls[,1], ls[,2])
ls = na.omit(data.frame(b1=progress[,col_k], b2=progress2[,col_k]))
cor(ls[,1], ls[,2])

grid.newpage() # 生成新画布
pushViewport(viewport(layout = grid.layout(2,3))) # 指定一个2乘3的layout
print(plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) # 画
print(plot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) # 画
print(plot3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3)) # 画
print(plot5, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) # 画
print(plot6, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)) # 画
print(plot7, vp = viewport(layout.pos.row = 2, layout.pos.col = 3)) # 画
popViewport() # 退出当前位置
############################################################



kalmanAllRates = function(t_time, tL, foursu_exons, total_exons, total_introns, splitnum=100, JW=0.0001){
  
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
    bt_dlm = (at_dlm-dropFirst(pre$s)[,4])/(dropFirst(pre$s)[,3]-dropFirst(pre$s)[,1])
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

checkdata <- data.frame(labexon=labexon['5034',],totexon=totexon['5034',],totintr=totintr['5034',],
                        labexon2=labexon2['5034',],totexon2=totexon2['5034',],totintr2=totintr2['5034',])

testres1=kalmanAllRates(t_time, tL, rbind(labexon['5034',]), rbind(totexon['5034',]), rbind(totintr['5034',]))
testres2=kalmanAllRates(t_time, tL, rbind(labexon2['5034',]), rbind(totexon2['5034',]), rbind2(totintr['5034',]))

# t4时刻 
trueres = c(-testres1$rates$ct[4],-testres1$rates$bt[4],testres1$rates$at[4])


x1=testres1$predata$gene[4,1]
x2=testres1$predata$gene[4,3] - x1
x1_dot=testres1$predata$gene[4,2]
x2_dot=testres1$predata$gene[4,4]

x1*(-testres1$rates$ct[4])+testres1$rates$at[4]
x2*(-testres1$rates$bt[4])+testres1$rates$at[4]


u=1
W = -matrix(c(x1**2, x1*x2, 0, 0, x1*u, 0,
                x2*x1, x2**2, 0, 0, x2*u, 0,
                0, 0, x1**2, x1*x2, 0, x1*u,
                0, 0, x2*x1, x2**2, 0, x2*u,
                u*x1, u*x2, 0, 0, u**2, 0,
                0, 0, u*x1, u*x2, 0, u**2), nrow=6)

I = matrix(c(x1*x1_dot, x2*x1_dot, x1*x2_dot, x2*x2_dot, u*x1_dot, u*x2_dot))

W = -matrix(c(x1**2,  0, x1*u, 
              0,  x2**2, x2*u, 
              u*x1, u*x2, 2*u**2), nrow=3, byrow = TRUE)

I = matrix(c(x1*x1_dot, x2*x2_dot, u*x1_dot+u*x2_dot))


#lamda = 5
#rou = 800 
#V =  rnorm(3)
#step = 0.00000001
#ts=c()
#counters = 1000
#parm = c()
#E = c()
#W%*%trueres + I
#W%*%V + I
#for(i in 1:counters){
#  dv = lamda*(rou**2-V**2)*(W%*%V + I)/(2*rou)
#  V = V + dv * step
#  parm = rbind(parm, t(V))
#  E = c(E, -0.5 * t(V)%*%W%*%V - t(V)%*%I)
#  ts = c(ts, i)
#}
#
#par(mfrow = c(1, 4))
#plot(ts, parm[,1])
#plot(ts, parm[,2])
#plot(ts, parm[,3])
#plot(ts, E)
#print(t(V))
#print(trueres)
#dv = get_dv(rou, lamda, V, W, I)
#V = V + dv.T * step






checkpoints = data.frame(na.omit(data.frame(b1=dlmbt[,col_k], b2=dlmbt2[,col_k])))
checkpoints = checkpoints[order(checkpoints[,1],checkpoints[,2],decreasing = TRUE),]
head(checkpoints)

checkdata <- data.frame(labexon=labexon['10',],totexon=totexon['10',],totintr=totintr['10',],
                        labexon2=labexon2['10',],totexon2=totexon2['10',],totintr2=totintr2['10',])

checkdata$id <- t_time
checkdata = melt(checkdata, id.vars = 'id')
ggplot(checkdata, aes(id, value)) + geom_line(aes(color=variable), size=1)



dlmrates <- dlmRates(t_time, tL, my_dlm, rbind(labexon['257396',]), rbind(totexon['257396',]), rbind(totintr['257396',]))
dlmrates2 <- dlmRates(t_time, tL, my_dlm, rbind(labexon2['257396',]), rbind(totexon2['257396',]), rbind(totintr2['257396',]))

JGG <- matrix(c(0,1,0,0),nrow=2,byrow=TRUE)
GG <- matrix(c(1,1,0,1),nrow=2,byrow=TRUE)
FF <- matrix(c(1,0), nrow=1)
m0 <- c(0,0)
C0 <- diag(1,2)
W <- diag(0.0001,2)
V <- 1
X = matrix(dert, nrow = num_data)
X = cbind(X, rep(0,num_data))
X = cbind(X, rep(1,num_data))
my_dlm <- dlm(X=X, FF=FF,V=V, JGG=JGG, GG=GG, W=W,m0=m0, C0=C0)

TT1=totexon['257396',]
PT1=totintr['257396',]
TL1=labexon['257396',]
pre_TT1 <- dlmSmooth(TT1,my_dlm)
pre_PT1 <- dlmSmooth(PT1,my_dlm)


dlmRes <- data.frame(dropFirst(pre_TT1$s), dropFirst(pre_PT1$s), TT1,PT1)
colnames(dlmRes) = c('pre_TT','pre_TT_dot', 'pre_PT', 'pre_PT_dot', 'TT', 'PT')
dlmRes$id <- t_time
dlmRes = melt(dlmRes, id.vars = 'id')
plot1 <- ggplot(dlmRes, aes(id, value)) + geom_line(aes(color=variable), size=1)

# 计算速率
at_dlm = TL1/tL
ct_dlm = (at_dlm-dropFirst(pre_PT1$s)[,2])/dropFirst(pre_PT1$s)[,1]
bt_dlm = (at_dlm-dropFirst(pre_TT1$s)[,2])/(dropFirst(pre_TT1$s)[,1]-dropFirst(pre_PT1$s)[,1])

dlmRates = data.frame(at_dlm,ct_dlm,bt_dlm)
dlmRates$id = t_time
dlmRates = melt(dlmRates, id.vars = 'id')
plot2 <- ggplot(dlmRates, aes(id, value)) + geom_line(aes(color=variable), size=1)

TT1=totexon2['257396',]
PT1=totintr2['257396',]
TL1=labexon2['257396',]
pre_TT1 <- dlmFilter(TT1,my_dlm)
pre_PT1 <- dlmFilter(PT1,my_dlm)


dlmRes <- data.frame(dropFirst(pre_TT1$m), dropFirst(pre_PT1$m), TT1,PT1)
colnames(dlmRes) = c('pre_TT','pre_TT_dot', 'pre_PT', 'pre_PT_dot', 'TT', 'PT')
dlmRes$id <- t_time
dlmRes = melt(dlmRes, id.vars = 'id')
plot3 <- ggplot(dlmRes, aes(id, value)) + geom_line(aes(color=variable), size=1)

# 计算速率
at_dlm = TL1/tL
ct_dlm = (at_dlm-dropFirst(pre_PT1$m)[,2])/dropFirst(pre_PT1$m)[,1]
bt_dlm = (at_dlm-dropFirst(pre_TT1$m)[,2])/(dropFirst(pre_TT1$m)[,1]-dropFirst(pre_PT1$m)[,1])

dlmRates = data.frame(at_dlm,ct_dlm,bt_dlm)
dlmRates$id = t_time
dlmRates = melt(dlmRates, id.vars = 'id')
plot4 <- ggplot(dlmRates, aes(id, value)) + geom_line(aes(color=variable), size=1)
##################################################

grid.newpage() # 生成新画布
pushViewport(viewport(layout = grid.layout(2,2))) # 指定一个2乘3的layout
print(plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) # 画
print(plot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) # 画
print(plot3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) # 画
print(plot4, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)) # 画
