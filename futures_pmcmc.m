%% run particle gibbs
% the following aims to sample from the posterior distribution over
% changepoint locations.
clear
rng(2)
addpath(genpath('util'))
addpath(genpath('csmc_cp'));
addpath(genpath('smc_cp'));
addpath(genpath('fearnhead_cp'));

% experiment settings
n_quadrature_points = 1e2;
constrain_rho = false;

experiment_dir = 'experiments/futures_pmcmc';
N = 10;
M = 1e4;

distribution_names = {'ghyp','sstd'};            
GDPPrior_arr = [false,true];


%% perform experiment on each of the data sources

for dist_indx=1:2
    
    % iterate through each of the distributions
    distribution_name = distribution_names{dist_indx};
    data_fname = sprintf('data/yahoo_data/largestocks2007%s.mat',distribution_name);
    load(data_fname);
    X = X_mat';
    [n,T] = size(X);

    for GDPPrior=GDPPrior_arr

        z0 = randn(1,T);
        pGeo0 = 0.5;
        tau = [0];

        % save settings and data to struct
        params.GDPPrior = GDPPrior;
        params.n_quadrature_points = n_quadrature_points;
        params.constrain_rho = constrain_rho;
        params.X = X;
        params.z = z0;
        params.pGeo = pGeo0;

        a_prior = 1;
        b_prior = 1;
        params.a_prior = a_prior;
        params.b_prior = b_prior;

        tau_collect = [];
        like_ests = zeros(1,M);
        pGeo_collect = zeros(1,M);
        corr_mat_collect =  zeros(n,n,M);
        z_collect = zeros(M,T);
        tic;
        for m=1:M
            disp(m)

            % update all params
            [tau,params,corr_mat,joint_loglik,rho_record] = particleGibbs(N,tau,params);

            % save results
            tau_collect = [tau_collect,tau];
            pGeo_collect(m) = params.pGeo;
            z_collect(m,:) = params.z;
            corr_mat_collect(:,:,m) = corr_mat;
        end
        time_taken = toc;
        
        save_fname = sprintf('%s/distribution_%s_gdpprior_%s.mat',...
            experiment_dir,distribution_name,string(GDPPrior));
        save(save_fname);

    end
end









