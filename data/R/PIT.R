require('rugarch')



logprices = readMat('logprices.mat')$logprices
n = dim(logprices)[2]
Tfinal = dim(logprices)[1]

u = matrix(0,Tfinal,n)

spec = ugarchspec()
for (i in (1:n)[-c(14,24)]){

    fit = ugarchfit(spec = spec, data = logprices[,i],solver = "nlminb")

  
  u[,i] = pit(fit)
}


writeMat('ulogprices.mat',u=u)