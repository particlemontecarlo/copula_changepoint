rm(list = ls())
require('rugarch')
library(R.matlab)
if (!require("quantmod")) {
  install.packages("quantmod")
  library(quantmod)
}

start <- as.Date("2007-01-01")
end <- as.Date("2009-01-01")


tickers <- c('XOM','MSFT','GE','IBM','CVX','JNJ','PG','WMT','BAC','GOOG')
n_tickers <- length(tickers)

getSymbols(tickers, src = "yahoo", from = start, to = end)

series_lst <- list(XOM,MSFT,GE,IBM,CVX,JNJ,PG,WMT,BAC,GOOG)
index_test <- index(series_lst[[1]])
n_dates <- length(index_test)
res <- rep(NA,n_tickers)
for(i in 2:n_tickers){
  res[i] <- all(index_test==index(series_lst[[i]]))
}



# get the PITs for each of the time series
distribution_models <- c('sstd','ghyp')
shapResults <- matrix(NA,length(distribution_models),n_tickers)
for(i in 1:length(distribution_models)){
  distribution_model <- distribution_models[i]
  ugspec = ugarchspec(distribution.model = distribution_model)
  
  shapResult = rep(NA,1,n_tickers)
  X_mat <- matrix(NA,n_dates-1,n_tickers)
  for(j in 1:n_tickers){
    series <- series_lst[[j]]
    names(series) <- c('Open','High','Low','Close','Volume','Adjusted')
    prices <- series$Close
    
    logp = log(prices)
    logreturns = diff(logp)[2:length(prices)]
    
    fit = ugarchfit(spec = ugspec, data = logreturns,solver = "solnp")
    #print(fit)
    u = pit(fit)
    
    Xnew = qnorm(u)
    shapResult[j] = shapiro.test(as.numeric(Xnew))$p.value
    X_mat[,j] <- Xnew
  }
  shapResults[i,] <- shapResult
  filename <- sprintf('largestocks2007%s.mat', distribution_model)
  writeMat(filename, X_mat=X_mat)
}
write.csv(index_test,file='dates.csv')
write.csv(tickers,file='tickers.csv')



