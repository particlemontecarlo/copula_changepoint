function [tau,params,corr_mat,joint_loglik,rho_record] = particleGibbs(N,tau_star,params)
%PARTICLEGIBBS Performs particle Gibbs for changepoint locations and
%generic Gibbs updates for the other parameters. Optionally add 
%(...,'ConstrainRho',true) to the input arguments so that inference is
%performed without label switching but with an identifiability constraint
%on the first component

% unpack parameters
X               = params.X;
z               = params.z;
a_prior         = params.a_prior;
b_prior         = params.b_prior;
GDPPrior        = params.GDPPrior;
constrain_rho   = params.constrain_rho;

[n,T] = size(params.X);

% cSMC to sample kernel leaving the full conditional over changepoints
% invariant
[SS_all,log_W_all,log_W_bar_all] = forwardFilteringCSMC(N,tau_star,params);
[ tau ] = bwdsSampling(params,SS_all,log_W_all)
[like_est] = get_likelihood_est(log_W_bar_all);

% sample factors and factor loadings conditional on changepoint
% locations
K = length(tau);
rhos = zeros(n,K);
t_vals = zeros(1,T);
rho_record = zeros(1,T);
for k=1:K
    % get the segment ends
    tstart = tau(k)+1;
    if k==K,tend = T;
    else, tend = tau(k+1); end

    % sample the rhos approximately
    for i = 1:n
        X_vals = X(i,tstart:tend);
        z_vals = z(tstart:tend);
        
        % optionally constrain rho on the first coordinate
        if constrain_rho && (i==1)
            rho_grid = linspace(0,1-1e-6,1e3);
        else
            rho_grid = linspace(-1+1e-6,1-1e-6,1e3);
        end
        
        log_prho = logpXiRhoiZ(X_vals,z_vals,rho_grid,GDPPrior);
        rho_cdf = cumsum(exp(log_prho));
        u_val = rand(1)*rho_cdf(end);
        rhos(i,k) = rho_grid(find(rho_cdf>=u_val,1));
    end
    
    if constrain_rho, assert(rhos(1,k)>=0), end

    % sample z, flipping within a segment if not constraining rho
    rhoj = rhos(:,k);
    varZpost = (sum(rhoj.^2./(1-rhoj.^2)) + 1 )^-1;
    
    if ~constrain_rho, s = -1+2*binornd(1,0.5,1,1);
    else, s = 1;
    end
    
    for t=tstart:tend
        muZpostt =  sum(rhoj.*X(:,t)./(1-rhoj.^2));
        z_samp = normrnd(muZpostt*varZpost,varZpost^0.5,1);
        z(t) = s*z_samp;
        rho_record(t) = rhoj(2);
        if t==1
            corr_mat = rhoj*rhoj';
        end
    end

    t_vals(tstart:tend) = tstart:tend;
end    
assert(all(t_vals==1:T));

% update params with new z value
params.z = z;

% update changepoint parameters conditional on changepoints
n_success = length(tau)-1;
params.pGeo = betarnd(a_prior + n_success,b_prior + T - n_success);

% log likelihood after marginalising out lambdas
joint_loglik = like_est + -0.5*(T*log(2*pi) + sum(z.^2)) - betalike([length(tau),1+T],params.pGeo);


end

