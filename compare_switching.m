%% run particle gibbs
% the following aims to sample from the posterior distribution over
% changepoint locations.
clear
rng(3)
addpath(genpath('util'))
addpath(genpath('csmc_cp'));
addpath(genpath('smc_cp'));
addpath(genpath('fearnhead_cp'));

% generate data
n = 10;
J = 2;
pGeo = 0.02;
segment_same_length = true;
GDPPrior = true;
[ X,K,rhosTrue,zTrue ] = generateCopulaData( n,J,pGeo,segment_same_length,GDPPrior );
corr_mat_true = rhosTrue(:,1)*rhosTrue(:,1)';


tau0 = cumsum(K);
[~,T] = size(zTrue);
n_quadrature_points = 1e2;
z = randn(1,T);
constrain_rho = true;


% save data to struct
params.X = X;
params.z = z;
params.GDPPrior = GDPPrior;
params.pGeo = pGeo;
params.n_quadrature_points = n_quadrature_points;



%%
n_vals = [10];
experiment_dir = 'experiments/compare_switching';
N = 10;
M = 1e5;



a_prior = 1;
b_prior = 1;
tau = 0:10:T;
params.z = randn(1,T);
params.a_prior = a_prior;
params.b_prior = b_prior;
params.constrain_rho = constrain_rho;

tau_collect = [tau];
like_ests = zeros(1,M);
pGeo_collect = zeros(1,M);
corr_mat_collect =  zeros(n,n,M);
z_collect = zeros(M,T);

tic;
for m=1:M
    disp(m)

    % update all params
    [tau,params,corr_mat] = particleGibbs(N,tau,params);

    % save results
    %like_ests(m) = like_est;
    tau_collect = [tau_collect,tau];
    pGeo_collect(m) = params.pGeo;
    z_collect(m,:) = params.z;
    corr_mat_collect(:,:,m) = corr_mat;
end
time_taken = toc;

save(sprintf('%s/compare_switching_n%i',experiment_dir,n))















