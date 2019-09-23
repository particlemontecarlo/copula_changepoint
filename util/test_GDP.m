%% script to test GDP sampling and density
% note the different parameterisations of the GDP distribution. We note
% that the default parameter values of (alpha,beta)=(3,1) in Murray 2011
% correspond to (alpha,xi)=(alpha,beta/alpha)=(3,1/3);

alpha = 3;
xi = 1/3;
n = 1e4;
X = sampleGDP(alpha,xi,n);

X_plt = linspace(-4,4);
[ p_plt ] = GDPdensity( X_plt,alpha,xi );

figure(1)
histogram(X,linspace(-4,4),'Normalization','probability');
hold all,plot(X_plt,0.095*p_plt),hold off