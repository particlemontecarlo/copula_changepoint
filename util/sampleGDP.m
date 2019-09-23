function [ res ] = sampleGDP( alpha,xi,n )
%SAMPLEGDP Samples generalised double Pareto distribution using Gibbs


pd = makedist('GeneralizedPareto','theta',0,'k',1/alpha,'sigma',xi);
res = random(pd,1,n);
random_sign = 2*(randi(2,1,n)-1)-1;
res = res.*random_sign;

end

