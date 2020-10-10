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

dirpath = 'E:\\wx\\2018下课题\\data'
rpkms_rep1=list()
rpkms_rep1$foursu_exons = read.table(paste(dirpath, '\\expression_gene.4sU.M.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$foursu_introns = read.table(paste(dirpath, '\\expression_gene.4sU.P.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$total_exons = read.table(paste(dirpath, '\\expression_gene.total.M.txt',sep=''),sep='\t', header = TRUE)
rpkms_rep1$total_introns = read.table(paste(dirpath, '\\expression_gene.total.P.txt',sep=''),sep='\t', header = TRUE)


rownames(rpkms_rep1$foursu_exons) = rpkms_rep1$foursu_exons[,1]
rownames(rpkms_rep1$foursu_introns) = rpkms_rep1$foursu_introns[,1]
rownames(rpkms_rep1$total_exons) = rpkms_rep1$total_exons[,1]
rownames(rpkms_rep1$total_introns) = rpkms_rep1$total_introns[,1]

rpkms_rep1$foursu_exons = rpkms_rep1$foursu_exons[,-1]
rpkms_rep1$foursu_introns = rpkms_rep1$foursu_introns[,-1]
rpkms_rep1$total_exons = rpkms_rep1$total_exons[,-1]
rpkms_rep1$total_introns = rpkms_rep1$total_introns[,-1]

# gene2ensimble
gene2ensimble = read.table(paste(dirpath,'\\id2geneid.txt',sep=''),sep='\t', header = TRUE)
gene2ensimble = as.matrix(gene2ensimble)[,1:2]

library("AnnotationDbi")
library("org.Hs.eg.db")


GENEID <- c(1,2,3,9,10)

nmid <- rownames(rpkms_rep1$foursu_exons)[1:5]

df <- data.frame(ENSEMBL,GENEID)
keytypes(org.Hs.eg.db)
# ENSG转换Symbol
df$symbol <- mapIds(org.Hs.eg.db,
                    keys=ENSEMBL,
                    column="ENSEMBL",
                    keytype="REFSEQ",
                    multiVals="first")





#px = read.table(paste(dirpath, '\\rates_gene.tx.txt',sep=''),sep='\t', header = TRUE)[,-1]
#tx = read.table(paste(dirpath, '\\rates_gene.px.txt',sep=''),sep='\t', header = TRUE)[,-1]
#dg = read.table(paste(dirpath, '\\rates_gene.dg.txt',sep=''),sep='\t', header = TRUE)[,-1]
############提取部分基因####################
labexon = rpkms_rep1$foursu_exons
labintr = rpkms_rep1$foursu_introns
totexon = rpkms_rep1$total_exons
totintr = rpkms_rep1$total_introns

genelist = intersect(rownames(labexon),rownames(totexon))
genelist = intersect(genelist,rownames(totintr))
genelist = intersect(genelist,rownames(labintr))

labexon = labexon[genelist,]
labintr = labintr[genelist,]
totexon = totexon[genelist,]
totintr = totintr[genelist,]

labexon = labexon[order(rownames(labexon),decreasing = TRUE),]
labintr = labintr[order(rownames(labintr),decreasing = TRUE),]
totexon = totexon[order(rownames(totexon),decreasing = TRUE),]
totintr = totintr[order(rownames(totintr),decreasing = TRUE),]



############INSPECT速率###################
t_time <- c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180)
tL <- 10
num_data = length(t_time)
dert = c(0, diff(t_time))

myrates <- newINSPEcT(t_time, tL, labexon, totexon,
                      labintr, totintr, BPPARAM=SerialParam())

##################################
# filter genes
ins_gene = rownames(ratesFirstGuess(myrates, 'synthesis'))
labexon_inspect = labexon[ins_gene,]
totexon_inspect = totexon[ins_gene,]
totintr_inspect = totintr[ins_gene,]
###############################################################################
# sort
labexon_inspect = labexon_inspect[order(labexon_inspect[,1],labexon_inspect[,2],
                                        labexon_inspect[,3],labexon_inspect[,4],labexon_inspect[,5],decreasing = TRUE),]

genename = rownames(labexon_inspect)
totexon_inspect = totexon_inspect[genename,]
totintr_inspect = totintr_inspect[genename,]

labexon_inspect = as.matrix(labexon_inspect)
totexon_inspect = as.matrix(totexon_inspect)
totintr_inspect = as.matrix(totintr_inspect)

res1 = estimateRates(t_time, tL, labexon_inspect, totexon_inspect+totintr_inspect, totintr_inspect, splitnum=100, JW=0.01)

###################################
filter.outers <- function(data, fold=1.5){
  filter_na = na.omit(data)
  eval_var=apply(filter_na, 1, var)
  fenwei = quantile(eval_var)
  # 去除异常值：将个体是距离底线或顶线的距离超过fold倍的箱体高度视为离群值。
  filter_var <- filter_na[(eval_var < (fenwei[4] + (fenwei[4]-fenwei[2])*fold) &
                             eval_var > (fenwei[2] - (fenwei[4]-fenwei[2])*fold)),]
  return(filter_var)
}
# dlm
dlmat = res1$at
dlmct = res1$ct
dlmbt = res1$bt

dlmat = filter.outers(dlmat)
dlmct = filter.outers(dlmct)
dlmbt = filter.outers(dlmbt)

synthe = ratesFirstGuess(myrates, 'synthesis')
progress = ratesFirstGuess(myrates, 'processing')
degradation = ratesFirstGuess(myrates, 'degradation')

gl = intersect(rownames(dlmat),rownames(dlmct))
gl = intersect(gl,rownames(dlmbt))
synthe = synthe[gl,]
progress = progress[gl,]
degradation = degradation[gl,]

dlmat = dlmat[gl, ]
dlmct = dlmct[gl, ]
dlmbt = dlmbt[gl, ]

################################

col_k = 10
plot1 <- ggplot(na.omit(data.frame(a1=dlmat[,col_k], a2=synthe[,col_k])), aes(a1, a2)) + geom_point()
plot2 <- ggplot(na.omit(data.frame(c1=dlmct[,col_k], c2=progress[,col_k])), aes(c1, c2)) + geom_point()
plot3 <- ggplot(na.omit(data.frame(b1=dlmbt[,col_k], b2=degradation[,col_k])), aes(b1, b2)) + geom_point() #+xlim(0,125)+ylim(0,100)

#plot5 <- ggplot(na.omit(data.frame(a1=synthe_sort[,col_k], a2=tx_drill_sort[,col_k])), aes(a1, a2)) + geom_point()
#plot6 <- ggplot(na.omit(data.frame(c1=progress_sort[,col_k], c2=px_drill_sort[,col_k])), aes(c1, c2)) + geom_point()
#plot7 <- ggplot(na.omit(data.frame(b1=degradation_sort[,col_k], b2=dg_drill_sort[,col_k])), aes(b1, b2)) + geom_point()

grid.newpage() # 生成新画布
pushViewport(viewport(layout = grid.layout(1,3))) # 指定一个2乘3的layout
print(plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) # 画
print(plot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2)) # 画
print(plot3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3)) # 画
#print(plot5, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) # 画
#print(plot6, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)) # 画
#print(plot7, vp = viewport(layout.pos.row = 2, layout.pos.col = 3)) # 画
popViewport() # 退出当前位置

##########################################
boxplot(dlmat)
#genelist = cbind(dlmat,apply(dlmat, 1, var))
#genelist = genelist[order(genelist[,14], decreasing = TRUE),]
#
#plot(genelist[1,], type='l', ylim=c((min(genelist)-1), (max(genelist)+1)))
#lines(genelist[2,], col='red')
#lines(genelist[3,], col='green')
#nameofgene = rownames(genelist)[c(1,2,3)]
lsr = c()
for(i in 1:10){
  profile_estimate = signalProfile(dlmat, dlmbt, dlmct,  ttime=t_time, topnum=2, loopnum=10, plotpic=TRUE)
  lsr = rbind(lsr, profile_estimate[,1])
}
profile_estimate = signalProfile2(dlmat,ttime=t_time, topnum=2, loopnum=10, plotpic=TRUE)

lsres = as.data.frame(t(lsr))
lsres$x = ttime
md = melt(lsres, id=c("x"))
ggplot(md, aes(x=x,y=value, color=variable))+geom_line()

lsr2 = c()
for(i in 1:10){
  profile_estimate2 = gensaProfile(dlmat, dlmbt, dlmct,  ttime=t_time, topnum=1, loopnum=10, plotpic=TRUE)
  lsr2 = rbind(lsr2, profile_estimate2)
}

lsres2 = as.data.frame(t(lsr2))
colnames(lsres2) = c(1,2,3)
lsres2$x = ttime
md2 = melt(lsres2, id=c("x"))
ggplot(md2, aes(x=x,y=value, color=variable))+geom_line()



ggplot(data.frame(y1=profile_estimate[,1],x=ttime, y2=profile_estimate2))+
  geom_smooth(aes(x=x,y=y2),method='loess')+
  geom_smooth(aes(x=x,y=y1),method='loess')
######################################################
# 估计gene
fai_t=profile_estimate[,1]
steadyrates = steadyRates(dlmat[1:20,], dlmbt[1:20,], dlmct[1:20,], t_time, profile_estimate[,1], loopnum=10, plotpic=FALSE)
num=19
plot(c(c(dlmat[num,]), dlmbt[num,], dlmct[num,]),type = 'o')
prea = steadyrates$a[num,1] + steadyrates$a[num,2]* fai_t * utFunc(ttime, steadyrates$a[num,3])
preb = steadyrates$b[num,1] + steadyrates$b[num,2]* fai_t * utFunc(ttime, steadyrates$b[num,3])
prec = steadyrates$c[num,1] + steadyrates$c[num,2]* fai_t * utFunc(ttime, steadyrates$c[num,3])
lines(c(prea,preb,prec), col='red',type = 'o')
plot(dlmct[num,],type = 'o')
lines(prec, col='red',type = 'o')


gensarates = gensaRates(dlmat[1:20,], dlmbt[1:20,], dlmct[1:20,], t_time, profile_estimate2, loopnum=10, plotpic=FALSE)
num=19
plot(c(c(dlmat[num,]), dlmbt[num,], dlmct[num,]),type = 'o')
prea = gensarates$a[num,1] + gensarates$a[num,2]* profile_estimate2 * utFunc(ttime, gensarates$a[num,3])
preb = gensarates$b[num,1] + gensarates$b[num,2]* profile_estimate2 * utFunc(ttime, gensarates$b[num,3])
prec = gensarates$c[num,1] + gensarates$c[num,2]* profile_estimate2 * utFunc(ttime, gensarates$c[num,3])
lines(c(prea,preb,prec), col='green',type = 'o')
plot(dlmct[num,],type = 'o')
lines(prec, col='green',type = 'o')

plot(gensarates$b[,1], steadyrates$b[,1])
##########################################################################################

signalProfile22 <-  function(rates_a, rates_b, rates_c, ttime, topnum=3, DriverGene=NULL ,loopnum=100, plotpic=FALSE){
  # 速率按照基因名字取交集
  rates_a = na.omit(as.matrix(rates_a))
  rates_b = na.omit(as.matrix(rates_b))
  rates_c = na.omit(as.matrix(rates_c))

  len = length(ttime)

  uniquegene = intersect(rownames(rates_a), rownames(rates_b))
  uniquegene = intersect(uniquegene, rownames(rates_c))
  rates_a = rates_a[uniquegene,]

  genelist_up = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, 1:len)))
  genelist_up = genelist_up[order(genelist_up[,(length(ttime)+1)], decreasing = TRUE),]

  genelist_down = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, len:1)))
  genelist_down = genelist_down[order(genelist_down[,(length(ttime)+1)], decreasing = TRUE),]

  genelist_up_down = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, c(1:(len%/%2),ceiling(len/2):1))))
  genelist_up_down = genelist_up_down[order(genelist_up_down[,(length(ttime)+1)], decreasing = TRUE),]

  genelist_down_up = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, c(ceiling(len/2):1, 1:(len%/%2)))))
  genelist_down_up = genelist_down_up[order(genelist_down_up[,(length(ttime)+1)], decreasing = TRUE),]

  if(is.null(DriverGene)){
    nameofgene = c(rownames(genelist_up)[1], rownames(genelist_down)[1],
                   rownames(genelist_up_down)[1],rownames(genelist_down_up)[1])
  }else{
    nameofgene = DriverGene
  }


  gene_rates = cbind(rates_a[nameofgene,], rates_b[nameofgene,])
  gene_rates = cbind(gene_rates, rates_c[nameofgene,])

  ttime <<- ttime
  gene_rates_data <<- gene_rates

  # 迭代loopnum次，计算fai
  resfai = c()
  #pb <- progress_bar$new(total = loopnum)
  pb <- progress_bar$new(
    format = "  running [:bar] :percent 执行时间 :elapsed 剩余: :eta ",
    total = loopnum, clear = FALSE, width= 120)
  for(i in 1:loopnum){
    lower = matrix(rep(c(0,0,0,-1,-1,-1,0,0,0, rep(0,length(ttime))), dim(gene_rates_data)[1]), nrow=dim(gene_rates_data)[1], byrow=TRUE)
    upper = matrix(rep(c(Inf,Inf,Inf,1,1,1,0.1,0.1,0.1, c(0.01,rep(Inf,(length(ttime)-1)))), dim(gene_rates_data)[1]), nrow=dim(gene_rates_data)[1], byrow=TRUE)
    newpar = setParm(rownames(gene_rates_data), ttime)
    tmpres = optim(par = c(newpar, newParameters(rep(1,(length(ttime))), lower = 0)),
                   energyFunc, NULL,
                   method = "L-BFGS-B",
                   lower = c(lower[,1:9], lower[1,10:(9+length(ttime))]),
                   upper = c(upper[,1:9], upper[1,10:(9+length(ttime))]))
    resfai = rbind(resfai,tmpres$par[(length(tmpres$par)-length(ttime)+1):length(tmpres$par)])

    pb$tick()
    Sys.sleep(1 / loopnum)
  }
  # 计算fai置信区间
  conffai = apply(resfai[,2:length(ttime)], 2, t.test)
  faidf = c()
  for(i in conffai){faidf = rbind(faidf,c(i$estimate,i$conf.int))}
  faidf = rbind(c(0,0,0), faidf)
  colnames(faidf) = c('mean', 'downlevel', 'uplevel')
  faidf = as.data.frame(faidf)
  faidf$x = ttime
  if(plotpic){
    g=ggplot(faidf)+
      geom_smooth(aes(x=x,y=mean),method='loess')
    print(g)
    #plot(ttime,faidf[,3], type = 'l',col='red')#, ylim=c(min(faidf), max(faidf)))
    #lines(ttime, faidf[,2], col='green')
    #lines(ttime, faidf[,1], col='brown')
  }

  return(faidf)
}


gensaProfile2 <-  function(rates_a, rates_b, rates_c, ttime, topnum=3, DriverGene=NULL ,loopnum=100, plotpic=FALSE){
  # 速率按照基因名字取交集
  rates_a = na.omit(as.matrix(rates_a))
  rates_b = na.omit(as.matrix(rates_b))
  rates_c = na.omit(as.matrix(rates_c))

  len = length(ttime)

  uniquegene = intersect(rownames(rates_a), rownames(rates_b))
  uniquegene = intersect(uniquegene, rownames(rates_c))
  rates_a = rates_a[uniquegene,]

  genelist_up = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, 1:len)))
  genelist_up = genelist_up[order(genelist_up[,(length(ttime)+1)], decreasing = TRUE),]

  genelist_down = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, len:1)))
  genelist_down = genelist_down[order(genelist_down[,(length(ttime)+1)], decreasing = TRUE),]

  genelist_up_down = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, c(1:(len%/%2),ceiling(len/2):1))))
  genelist_up_down = genelist_up_down[order(genelist_up_down[,(length(ttime)+1)], decreasing = TRUE),]

  genelist_down_up = cbind(rates_a,apply(rates_a, 1, function(x)cor(x, c(ceiling(len/2):1, 1:(len%/%2)))))
  genelist_down_up = genelist_down_up[order(genelist_down_up[,(length(ttime)+1)], decreasing = TRUE),]

  if(is.null(DriverGene)){
    nameofgene = c(rownames(genelist_up)[1], rownames(genelist_down)[1],
                   rownames(genelist_up_down)[1],rownames(genelist_down_up)[1])
  }else{
    nameofgene = DriverGene
  }


  gene_rates = cbind(rates_a[nameofgene,], rates_b[nameofgene,])
  gene_rates = cbind(gene_rates, rates_c[nameofgene,])

  ttime <<- ttime
  gene_rates_data <<- gene_rates

  lower = matrix(rep(c(0,0,0,-1,-1,-1,0,0,0, rep(0,length(ttime))), dim(gene_rates_data)[1]), nrow=dim(gene_rates_data)[1], byrow=TRUE)
  upper = matrix(rep(c(Inf,Inf,Inf,1,1,1,0.1,0.1,0.1, c(0.01,rep(Inf,(length(ttime)-1)))), dim(gene_rates_data)[1]), nrow=dim(gene_rates_data)[1], byrow=TRUE)
  newpar = setParm(rownames(gene_rates_data), ttime)
  # 计算大致搜索空间
  message('计算大致搜索空间')
  resfai = c()
  for(i in 1:loopnum){
    tmpres = optim(par = c(newpar, newParameters(rep(1,(length(ttime))), lower = 0)),
                   energyFunc, NULL,
                   method = "L-BFGS-B",
                   lower = c(lower[,1:9], lower[1,10:(9+length(ttime))]),
                   upper = c(upper[,1:9], upper[1,10:(9+length(ttime))]))
    resfai = rbind(resfai,tmpres$par)
  }
  up = apply(as.data.frame(resfai), 2, max)
  lower = matrix(rep(c(0,0,0,-1,-1,-1,0,0,0, rep(0,length(ttime))), dim(gene_rates_data)[1]), nrow=dim(gene_rates_data)[1], byrow=TRUE)
  upper = matrix(up[1:(length(tmpres$par)-length(ttime))], nrow=dim(gene_rates_data)[1], byrow=TRUE)
  message('开始估计')
  tmpres = GenSA(fn = energyFunc,
                 lower = c(lower[,1:9], lower[1,10:(9+length(ttime))]),
                 upper = c(upper[,1:9], up[(length(tmpres$par)-length(ttime)+1):length(tmpres$par)])+10)
  resfai = tmpres$par[(length(tmpres$par)-length(ttime)+1):length(tmpres$par)]

  return(resfai)
}

gensaRates <-  function(rates_a, rates_b, rates_c, ttime, fai_t, loopnum=100, plotpic=FALSE){
  res = list(a=c(),b=c(),c=c())
  fai <<- fai_t
  # 速率按照基因名字取交集
  uniquegene = intersect(rownames(rates_a), rownames(rates_b))
  uniquegene = intersect(uniquegene, rownames(rates_c))
  gene_rates_data = cbind(rates_a[uniquegene,], rates_b[uniquegene,])
  gene_rates_data = cbind(gene_rates_data, rates_c[uniquegene,])
  # 过滤空值
  gene_rates_data <<- na.omit(gene_rates_data)
  dimgene = dim(gene_rates_data)[1]
  up = max(gene_rates_data)+5
  for(g in 1:dimgene){
    ##################
    single_gene_data <<- gene_rates_data[g,1: length(ttime)]
    atmp = GenSA(fn = singleGeneEnergy,
                 lower = c(0,-1,0),
                 upper = c(up,1,ttime[length(ttime)]))
    atmpc = atmp$par
    ###################
    single_gene_data <<- gene_rates_data[g,(length(ttime)+1) : (2*length(ttime))]
    btmp =  GenSA(fn = singleGeneEnergy,
                  lower = c(0,-1,0),
                  upper = c(up,1,ttime[length(ttime)]))
    btmpc = btmp$par
    ##################
    single_gene_data <<- gene_rates_data[g,(2*length(ttime)+1): (3*length(ttime))]
    ctmp =  GenSA(fn = singleGeneEnergy,
                  lower = c(0,-1,0),
                  upper = c(up,1,ttime[length(ttime)]))
    ctmpc = ctmp$par

    res$a=rbind(res$a,c(atmpc))
    res$b=rbind(res$b,c(btmpc))
    res$c=rbind(res$c,c(ctmpc))
  }
  return(res)
}











