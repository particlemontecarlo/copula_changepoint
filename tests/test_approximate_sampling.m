%% estimates posterior predictive
% the following script estimates the posterior predictive for observations
% arising in the factor copula model. The methodology is the start of 
clear
rng(2)
addpath(genpath('util'))
addpath(genpath('csmc_cp'));
addpath(genpath('smc_cp'));
addpath(genpath('fearnhead_cp'));

% generate data
n = 5;
J = 1;
pGeo = 0.1;
segment_same_length = true;
GDPPrior = true;
[ X,K,rhosTrue,zTrue ] = generateCopulaData( n,J,pGeo,segment_same_length,GDPPrior );
[~,T] = size(X);
z = zTrue;
corr_mat_true = rhosTrue(:,1)*rhosTrue(:,1)';


% rhos using array fun
n_rho = 1e3;
rho_grid = linspace(-1+1e-6,1-1e-6,n_rho);
i = 2;
prho = exp(logpXiRhoiZ(X(i,:),z,rho_grid,GDPPrior));
plot(rho_grid,prho/sum(prho))


% rhos using posterior and prior
p_rho = zeros(1,n_rho);
for s=1:n_rho
    rho = rho_grid(s);
    logp_arr = zeros(1,T);
    for t=1:T
        var = 1-rho^2;
        logp_arr(t) = -0.5*(log(2*pi*var) + (X(i,t)-rho*z(t))^2/var);
            
    end
    log_prior = log( GDPdensity( rho/sqrt(1-rho^2),3,1/3 )) + ...
            -1.5*log(1-rho^2);
    p_rho(s) = exp(sum(logp_arr)+log_prior);
end

hold all,plot(rho_grid,p_rho/sum(p_rho)),hold off




