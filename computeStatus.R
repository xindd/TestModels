library(progress)

newParameters <- function(par, lower=NULL, upper=NULL, sd=2){
  new.par = c()
  #在原参数附近产生新的参数
  for(p in as.vector(par)){
    s = rnorm(1,mean = p, sd = sd)
    if(is.null(lower) && !is.null(upper)){ # 没有最小值有最大值
      while(s>upper){
        s = rnorm(1,mean = p, sd = sd)
      }
    }else if(!is.null(lower) && is.null(upper)){ # 有最小值没有最大值
      while(s<lower){
        s = rnorm(1,mean = p, sd = sd)
      }
    }else if(!is.null(lower) && !is.null(upper)){ # 既有最小值又有最大值
      while(s<lower || s>upper){
        s = rnorm(1,mean = p, sd = sd)
      }
    }

    new.par = c(new.par, s)
  }
  return(new.par)
}
utFunc = function(t, tao){
  ut = t - tao
  ut[which(ut<0)]=0
  ut[which(ut>=0)]=1
  return(ut)
}
# 需要定义观测数据gene_rates_data，时间序列ttime
# 返回误差能量值
energyFunc = function(paramaters){
  observes = gene_rates_data
  params = matrix(paramaters[1:(dim(observes)[1]*9)],ncol = 9)
  fai = paramaters[(dim(observes)[1]*9+1):length(paramaters)]
  t = ttime
  # 计算基因能量
  etmp=0
  len = length(t)
  for(g in 1:dim(observes)[1]){
    par = params[g, ]
    observe = observes[g,]
    uta = utFunc(t, par[7])
    utb = utFunc(t, par[8])
    utc = utFunc(t, par[9])

    error_a = observe[1:len] - (par[1] + par[4] * fai *uta) #/ (1 + exp(-(t-par[7]))))
    error_b = observe[(len+1):(2*len)] - (par[2] + par[5] * fai * utb) #/ (1 + exp(-(t-par[8]))))
    error_c = observe[(2*len+1):(3*len)] - (par[3] + par[6] * fai * utc) #/ (1 + exp(-(t-par[9]))))
    etmp = etmp + sum(error_a**2 + error_b**2 + error_c**2)
  }
  #energyG <<- c(energyG, etmp)
  #atG <<- c(atG, paramaters[1])
  return(etmp)
}

# 设置参数
setParm <- function(genelist, ttime, rate=1, S=1, tao=1, tao_sd=2){
  initParm = c()
  for(g in as.vector(genelist)){
    newpam = c(newParameters(rep(rate,3), lower = 0), # a,b,c>0
               newParameters(rep(S,3), lower = -1, upper = 1), # Sa,Sb,Sc
               newParameters(rep(tao,3), lower = 0,upper = ttime[length(ttime)], sd=tao_sd)) # 0<ta,tb,tc<T[end]
    initParm = rbind(initParm, newpam)
  }
  rownames(initParm) = genelist
  colnames(initParm) = c('as', 'bs', 'cs', 'Sa', 'Sb', 'Sc', 'ta', 'tb', 'tc')
  return(initParm)
}
# 估计刺激函数谱
signalProfile <- function(rates_a, rates_b, rates_c, ttime, loopnum=100, plotpic=FALSE){
  # 速率按照基因名字取交集
  uniquegene = intersect(rownames(rates_a), rownames(rates_b))
  uniquegene = intersect(uniquegene, rownames(rates_c))
  gene_rates_data = cbind(rates_a[uniquegene,], rates_b[uniquegene,])
  gene_rates_data = cbind(gene_rates_data, rates_c[uniquegene,])
  # 过滤空值
  gene_rates_data <<- na.omit(gene_rates_data)
  if(length(gene_rates_data) <= (length(ttime)*3)){
    message('过滤后基因数量少于1')
    stop()
  }
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
  colnames(faidf) = c('mean', 'down-level', 'up-level')
  if(plotpic){
    plot(ttime,faidf[,1], type = 'l',col='red')
    lines(ttime, faidf[,2], col='green')
    lines(ttime, faidf[,3], col='brown')
  }
  return(faidf)
}
# 单基因能量函数
singleGeneEnergy <- function(paramaters){
  observer = single_gene_data
  ut = utFunc(ttime, paramaters[3])
  error = observer - (paramaters[1] + paramaters[2] * fai * ut) #/ (1 + exp(-(t-par[7]))))
  return(sum(error**2))
}
# 估计基因稳态速率
steadyRates <- function(rates_a, rates_b, rates_c, ttime, fai_t, loopnum=100, plotpic=FALSE){
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
  pb <- progress_bar$new(
    format = "  running [:bar] :percent 执行时间 :elapsed 剩余: :eta ",
    total = dimgene, clear = FALSE, width= 120)
  for(g in 1:dimgene){
    ##################
    single_gene_data <<- gene_rates_data[g,1: length(ttime)]
    atmpc = c()
    for(i in 1:loopnum){
      atmp = optim(par = c(newParameters(c(1,1,1), lower = 0)),
                   singleGeneEnergy, NULL,
                   method = "L-BFGS-B",
                   lower = c(0,-1,0),
                   upper = c(Inf,1,ttime[length(ttime)]))
      atmpc = rbind(atmpc, atmp$par)
    }
    confat = apply(atmpc, 2, t.test)
    atdf = c()
    for(i in confat){atdf = rbind(atdf,c(i$estimate,i$conf.int))}
    ###################
    single_gene_data <<- gene_rates_data[g,(length(ttime)+1) : (2*length(ttime))]
    btmpc = c()
    for(i in 1:loopnum){
      btmp = optim(par = c(newParameters(c(1,1,1), lower = 0)),
                   singleGeneEnergy, NULL,
                   method = "L-BFGS-B",
                   lower = c(0,-1,0),
                   upper = c(Inf,1,ttime[length(ttime)]))
      btmpc = rbind(btmpc, btmp$par)
    }
    confbt = apply(btmpc, 2, t.test)
    btdf = c()
    for(i in confbt){btdf = rbind(btdf,c(i$estimate,i$conf.int))}
    ##################
    single_gene_data <<- gene_rates_data[g,(2*length(ttime)+1): (3*length(ttime))]
    ctmpc = c()
    for(i in 1:loopnum){
      ctmp = optim(par = c(newParameters(c(1,1,1), lower = 0)),
                   singleGeneEnergy, NULL,
                   method = "L-BFGS-B",
                   lower = c(0,-1,0),
                   upper = c(Inf,1,ttime[length(ttime)]))
      ctmpc = rbind(ctmpc, ctmp$par)
    }

    confct = apply(ctmpc, 2, t.test)
    ctdf = c()
    for(i in confct){ctdf = rbind(ctdf,c(i$estimate,i$conf.int))}

    res$a=rbind(res$a,c(atdf))
    res$b=rbind(res$b,c(btdf))
    res$c=rbind(res$c,c(ctdf))

    pb$tick()
    Sys.sleep(1 / dimgene)
  }
  return(res)
}

# simulation
createSim = function(genenum, tt, f, tao=0, tao_sd=2){
  gl = c()
  ratea = c()
  rateb = c()
  ratec = c()
  for(p in 1:genenum){gl = c(gl, paste('g',p, sep=''))}
  sp = setParm(gl, tt, rate = 10, tao = 0,tao_sd = tao_sd)
  for(g in gl){
    gene = sp[g,]
    test_a = gene[1] + gene[4]*f*utFunc(tt, gene[7])
    test_b = gene[2] + gene[5]*f*utFunc(tt, gene[8])
    test_c = gene[3] + gene[6]*f*utFunc(tt, gene[9])

    ratea = rbind(ratea, test_a)
    rateb = rbind(rateb, test_b)
    ratec = rbind(ratec, test_c)
  }
  rownames(ratea) = gl
  rownames(rateb) = gl
  rownames(ratec) = gl
  return(list(param=sp, ratea=ratea, rateb=rateb, ratec=ratec))
}
###############################################
#######################################
ttime = c(0,1,2,3,4)
f = c(0,5,9,4,3)
simRes = createSim(3, ttime, f, tao_sd = 0.1)
testRes = createSim(20, ttime, f, tao_sd = 2)
# 估计fai谱
fia_conf = signalProfile(testRes$ratea[1:3,], testRes$rateb[1:3,], testRes$ratec[1:3,], ttime, loopnum=10, plotpic=TRUE)
lines(ttime, f, col='black')

plot(ttime,fia_conf[,1], type = 'l',col='red')
lines(ttime, fia_conf[,2], col='green')
lines(ttime, fia_conf[,3], col='brown')
# 估计gene
fai_t=fia_conf[,1]
steadyrates = steadyRates(testRes$ratea, testRes$rateb, testRes$ratec, ttime, fai_t, loopnum=100, plotpic=FALSE)
plot(c(c(t(testRes$param)),f),type = 'o')
resvector = c()
for(i in 1:20){
  resvector = c(resvector, c(steadyrates$a[1,1],steadyrates$b[1,1],steadyrates$c[1,1],
                             steadyrates$a[1,2],steadyrates$b[1,2],steadyrates$c[1,2],
                             steadyrates$a[1,3],steadyrates$b[1,3],steadyrates$c[1,3]))
}
lines(c(resvector,fai_t), col='red')


ut = ttime-steadyrates$a[1,3]
ut[which(ut<0)]=0
ut[which(ut>0)]=1


###############################################
