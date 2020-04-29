library(quantmod)
library(RConics)
library(PerformanceAnalytics)
library(moments)

BS <-
  function(S, K, T, r, sig, type="C"){
    d1 <- (log(S/K) + (r + sig^2/2)*T) / (sig*sqrt(T))
    d2 <- d1 - sig*sqrt(T)
    if(type=="C"){
      value <- S*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
    }
    if(type=="P"){
      value <- K*exp(-r*T)*pnorm(-d2) - S*pnorm(-d1)
    }
    return(value)
  }

logdiff = function(d){
  len = length(d[,1])
  rets =  d[,]
  for(j in 1:length(rets[1,])){
    rets[1:len,j] = diff(log(rets[,j]))
  }
  rets = rets[2:len,]
  return(rets)
}
OPToption = function(ra,price,mu,sigma,r_free){
  f = mu/(sigma*sigma)/ra
  
  k = c()
  for(i in 1:floor(price)){
    Bsmc = BS(price,i,1,r_free,sigma,type = "C")
    Intrinsic = price - i
    Extrinsic = Bsmc - Intrinsic
    Interest = Extrinsic/i
    #print(Extrinsic)
    Leverage = price/(price - i)
    #print((f - Interest/(sigma*sigma)) - Leverage)
    if(((f - Interest/(ra*sigma*sigma)) - Leverage) <= 0.02){
      
      k = rbind(i,Bsmc,Leverage)
      break
    }
  }
  return(k)
}

dataset = xts()
data <- getSymbols("^GSPC",from = "1989-01-01", to = "2019-11-02", src = 'yahoo')
dataset = Ad(get("GSPC"))
#data = getSymbols("^VIX",from = "1987-01-01", to = "2019-10-23", src = 'yahoo')
#dataset = merge(dataset,Ad(get("VIX")))
data <- getSymbols("AAA",from = "1989 - 01-01",to = "2019-11-02", src = "FRED")
dataset2 = xts()
dataset2 = (get("AAA"))
data <- getSymbols("DGS1",from = "1989 - 01-01",to = "2019-11-02", src = "FRED")
dataset2 = merge(dataset2,(get("DGS1")))
dataset2 = dataset2[(7654:length(dataset2[,1])),]
dataset2[1,2] = 9.11
class(dataset2)
class(dataset)
ddd = merge(dataset,dataset2)
ddd = ddd[complete.cases(ddd)]
for(i in 1:length(dataset2[,2])){
  if(is.na(dataset2[i,2]) == TRUE){
    dataset2[i,2] = dataset2[(i-1),2]
  }
  if(is.na(dataset2[i,1]) == TRUE){
    dataset2[i,1] = dataset2[(i-1),1]
  }
}
ddd = merge(dataset,dataset2)
ddd = ddd[complete.cases(ddd)]
 #initial value of bond portfolio.

#we want to trade only AAA rated bonds in the portfolio however returns are allowed to be traded
#in an attempt to maximize growth, how should we trade?
#we can either keep it in AAA rated bonds or we can trade in the SP500
#you get penalized for losing basically
#we need an algorithm to get best option to buy


Lrets = logdiff(ddd$GSPC.Adjusted)
ddd2 = merge(ddd,Lrets)
ddd2 = ddd2[complete.cases(ddd2)]


lever = c(0)
lever = as.matrix(lever)
ra = 2 #risk aversion parameter aka half kelly = 2 full kelly = 1

pVAR = rbind(1000000,1000000)
Pkelly = rbind(1000000,1000000)
Portfolio_value = c(1000000) #starting portfolio value
Portfolio_value = as.matrix(Portfolio_value)
rownames(Portfolio_value)[1] = c(as.character(index(ddd2[(252*2),2]))) #portfolio return dates we are denoting returns for that fiscal year, if year end starts at 1991 then the returns are for year 1990
rownames(lever)[1] = c(as.character(index(ddd2[(252*2),2])))
Portfolio_value = rbind(Portfolio_value,Portfolio_value[1]*(1+as.numeric(ddd2[1,2])/100)) #need 1 step to get profits
lever = rbind(lever,0)
rownames(lever)[2] = c(as.character(index(ddd2[(252*3),2])))
rownames(Portfolio_value)[2] = c(as.character(index(ddd2[(252*3),2])))
Benchmark = Portfolio_value
SPX = c(ddd2[(252*2),1],ddd2[(252*3),1])
rownames(SPX) = c(rownames(Benchmark))


checks = c()
for(i in 1:(floor(length(ddd2[,1])/252)- 3 )){
  ts = 1 + 252*(i-1)
  tf = 252*3 + 1 + 252*(i-1)
  
  
  
  m = mean(ddd2[(ts:tf),4])
  s = sd(ddd2[(ts:tf),4])
  
  var = 0.2 #value at risk factor, 0 < var < 1 denotes how much of portfolio you risk in market
  BSVAR = BS(as.numeric(ddd2[tf,1]),as.numeric(ddd2[tf,1])*(1-var),1,as.numeric(ddd2[tf,3])/100,s,type = "C")#Optionsprice where we keep VAR at 20% of portfolio
  payVAR = as.numeric(ddd2[(tf + 252),1]) - as.numeric(ddd2[tf,1])*(1-var)
  pVAR = rbind(pVAR,((1-var)*pVAR[i+1]*(1 + as.numeric(ddd2[tf,2])/100) + var*max(payVAR,0)*pVAR[i+1]/BSVAR))
  #pVAR = rbind(pVAR,((1-var)*pVAR[i+1]*(1 + as.numeric(ddd2[tf,3])/100) + var*max(payVAR,0)*pVAR[i+1]/BSVAR))
  #pVAR = rbind(pVAR,((1-var)*pVAR[i+1]*(1) + var*max(payVAR,0)*pVAR[i+1]/BSVAR))
  
  f = (m*252 - as.numeric(ddd2[tf,3])/100)/(s*s*252)
  play_money = 0
  if(Portfolio_value[i+1] - Portfolio_value[i] > 0){ #the trader can only trade last years profits, trader is punished by not being able to trade if returns are too low
    play_money = Portfolio_value[i+1] - Portfolio_value[i]
  }
  else{
    play_money = 0
  }
  ppp = c()
  if( f > 1){ 
    ppp = OPToption(ra,as.numeric(ddd2[tf,1]),m*252,s*sqrt(252),as.numeric(ddd2[tf,3])/100)
    f = ppp[3]
    payoff_option = max((ddd2[(tf + 252),1] - ppp[1])*1*(play_money/ppp[2]),0)
    pkelly_option = max((ddd2[(tf + 252),1] - ppp[1])*1*(Pkelly[i+1]/ppp[2]),0)
  }
  else if(f < 0){ #trader can only be long. i dont have a system for buying puts, so trader is long AAA bonds in this scenario
    f = 0
    play_money = 0
    payoff_option = 0
    pkelly_option = Pkelly[i+1]
  }
  else{ #if leverage is less than 1 trader can choose to just be long shares of SPX at some fraction of portfolio profits and be long AAA bonds with the rest
    f = f/ra
    play_money = f*play_money
    
    payoff_option = play_money + play_money*sum((ddd2[(tf:(tf + 252)) , 4]))
    pkelly_option = Pkelly[i+1] + Pkelly[i+1]*sum((ddd2[(tf:(tf + 252)) , 4]))
  }
 
 weight = c(1 - play_money/Portfolio_value[i + 1], play_money/Portfolio_value[i + 1])
 
 Pkelly = rbind(Pkelly,pkelly_option)
 total_return = weight[1]*Portfolio_value[i + 1]*(1+ddd2[tf,2]/100) + max(payoff_option,0)
 net = total_return - Portfolio_value[i+1] #portfolio returns
 net_option = payoff_option - play_money #cash returns on kelly bets
 netpercent = log(payoff_option/play_money) #log returns of kelly bets
 Portfolio_value = rbind(Portfolio_value, as.numeric(total_return)) #update portfolio value
 
 x = c(play_money,f,weight,payoff_option,net,net_option,netpercent)
 x = t(x)
 checks = rbind(checks,x)
 lever = rbind(lever,f)
 
 Benchmark = rbind(Benchmark, as.numeric(Benchmark[i+1]*(1+ddd2[tf,2]/100)))
 SPX = rbind(SPX,ddd2[tf + 252,1])
 rownames(SPX)[i+2] = c(as.character(index(ddd2[(tf + 252),2])))
 rownames(Benchmark)[i + 2] = c(as.character(index(ddd2[(tf + 252),2])))
 rownames(Portfolio_value)[i + 2] = c(as.character(index(ddd2[(tf + 252),2])))
 rownames(lever)[i + 2] = c(as.character(index(ddd2[(tf + 252),2])))
}
colnames(checks) = c("play_money","Kelly_factor/2","weight_AAA","Weight_Kelly","payoff_option","net_returns","Kelly_returns","Kelly_returns_%")
#plot(cumsum(Lrets),type = "l")
Pkelly = as.matrix(Pkelly)
rownames(Pkelly) = c(rownames(Portfolio_value))
rownames(pVAR) = c(rownames(Portfolio_value))
lll = log(as.xts(as.matrix(Portfolio_value)))
lll = diff(lll)
lll[1] = 0
bbb = log(as.xts(as.matrix(Benchmark)))
bbb = diff(bbb)
bbb[1] = 0
ccc = log(as.xts(as.matrix(SPX)))
ccc = diff(ccc)
ccc[1] = 0
ggg = diff(log(as.xts(Pkelly)))
ggg[1] = 0
qqq = diff(log(as.xts(pVAR)))
qqq[1] = 0

plot(as.xts(lever), type = "l")
plot(cumsum(merge(ccc,bbb,lll,ggg)),type = "l")
x = c(sd(ccc)/sd(lll),sd(lll - bbb))

plot(cumsum(merge(qqq[3:length(qqq)],ccc[3:length(ccc)])))
qqq[3:length(qqq)] - ccc[3:length(ccc)]