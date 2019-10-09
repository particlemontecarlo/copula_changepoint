%% run particle gibbs
% the following aims to sample from the posterior distribution over
% changepoint locations.
clear
rng(1)
addpath(genpath('util'))
addpath(genpath('csmc_cp'));
addpath(genpath('smc_cp'));
addpath(genpath('fearnhead_cp'));

% generate data
n = 15;
J = 5;
pGeo = 0.05;
segment_same_length = true;
GDPPrior = false;
[ X,K,rhosTrue,zTrue ] = generateCopulaData( n,J,pGeo,segment_same_length,GDPPrior );
corr_mat_true = rhosTrue(:,1)*rhosTrue(:,1)';

tau0 = cumsum(K);
[~,T] = size(zTrue);
n_quadrature_points = 1e2;
z = zTrue;

% save data to struct
params.X = X;
params.z = randn(1,T);
params.GDPPrior = GDPPrior;
params.pGeo = pGeo;
params.n_quadrature_points = n_quadrature_points;
params.constrain_rho = false;


%%
n_vals = [5,10,15];
experiment_dir = 'experiments/posterior_concentration';
if ~exist(experiment_dir, 'dir'),mkdir(experiment_dir),end

N = 10;
M = 30000;


a_prior = 1;
b_prior = 1;
tau = [0];

params.a_prior = a_prior;
params.b_prior = b_prior;

for ss=1:length(n_vals)
    n = n_vals(ss);
    params.X = X(1:n,:);

    tau_collect = [];
    like_ests = zeros(1,M);
    pGeo_collect = zeros(1,M);
    corr_mat_collect =  zeros(n,n,M);
    z_collect = zeros(M,T);
    rho_record_collect = zeros(M,T);

    tic;
    for m=1:M
        disp(m)

        % update all params
        [tau,params,corr_mat,joint_loglik,rho_record] = particleGibbs(N,tau,params);

        % save results
        %like_ests(m) = like_est;
        tau_collect = [tau_collect,tau];
        pGeo_collect(m) = params.pGeo;
        z_collect(m,:) = params.z;
        corr_mat_collect(:,:,m) = corr_mat;
        rho_record_collect(m,:) = rho_record;
    end
    time_taken = toc;
    
    save(sprintf('%s/posterior_concentration_n%i',experiment_dir,n))
end














