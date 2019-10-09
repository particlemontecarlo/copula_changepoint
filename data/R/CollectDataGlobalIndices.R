# require R timeseries library
require('zoo')
require('rugarch')
require('R.matlab')

#tickNames = c('BO','CC','CT','DA','GC','HG','JO','KC','KW','LB')
#tickNamesOriginal = c('ZZ','W_','HG','KW','NG','HO','LC')
#tickNames = tickNamesOriginal

#Passed the test:
# HG (Copper),KW (Wheat), NG (Natural gas)

tickNames10working = c('ZZ','W_','HG','KW','NG','FA','LC','GC','JO','KC') 

badList = c('AP','BC','BG','MW','SM','GI','LB')
tickNames20Additional1 = c('AD','AX','BN','CA','CB','CL','CR','DJ','DT','DX','EC','EN','ER','ES','FA','FC','FN','FX','GS')
tickNames20Additional2 = c('HS','JN','LH','LX','MD','MP','ND','NK','NR','PA','PL','RB','RL','SB','SC','SF','SI','ZB')
tickNamesAdditional = c('SN','SP','SS','TA','TD','UZ','XU','XX','YM')


tickNamesGlobalIndices = c('DJI','FTSE','GSPC','NDX','SSMI')


tickNames = tickNamesGlobalIndices#tickNamesGlobalIndices#c(tickNames10working,tickNames20Additional1,tickNames20Additional2,tickNamesAdditional)
N = length(tickNames)

ugspec = ugarchspec(distribution.model = "std")


urlData1 = "http://hapi.oxford-man.ox.ac.uk/hapi37.php?t=y&s=^"
urlData2 = "&f=close&d1=2005-01-01&d2=2010-01-01" 

Xall = zoo();


urlData = paste(urlData1,tickNames[1],urlData2,sep='')
prices <- read.zoo(urlData, header = TRUE, sep = ",", format = "%Y/%m/%d")
pricesAll = prices/as.numeric(prices[1])

logreturns = diff(log(prices))

fit = ugarchfit(spec = ugspec, data = logreturns,solver = "solnp")

u = pit(fit)
X = qnorm(u)

shapResult = rep(NA,1,N)
shapResult[1] = shapiro.test(as.numeric(X))$p.value
for (i in 2:N){

  urlData = paste(urlData1,tickNames[i],urlData2,sep='')
  prices <- read.zoo(urlData, header = TRUE, sep = ",", format = "%Y/%m/%d")
  pricesAll = merge(pricesAll,prices/as.numeric(prices[1]))
  
  logp = log(prices)
  logreturns = diff(logp)
  
  fit = ugarchfit(spec = ugspec, data = logreturns,solver = "solnp")
  #print(fit)
  u = pit(fit)
  
  Xnew = qnorm(u)
  shapResult[i] = shapiro.test(as.numeric(Xnew))$p.value
  X = merge(X,Xnew,all=c(T,T))
  print(i/N)
}

passed = shapResult>0.05
highNAs = colSums(is.na(X))>50
#passed = passed & !highNAs
passed = rep(T,1,N)

shapResultPost = shapResult[passed]
tickNamesPost = tickNames[passed]
XPost = X[,passed]

NAs = which(is.na(XPost),arr.ind=T)
rowNAs = unique(NAs[,1])
XPost = XPost[-rowNAs,]

Nnew = dim(XPost)[2]
Y = as.data.frame(XPost)
writeMat(paste('DataReadyStdMarginal012006to012010',Nnew,'series.mat'),X=Y,dates = rownames(Y),Names = tickNames,shapResult=shapResult)




