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
n = 10;
J = 3;
pGeo = 1/20;
segment_same_length = true;
GDPPrior = false;
fraction_change = 1;
[ X1,K,rhosTrue,zTrue ] = generateCopulaDataLowDimCP( n,J,pGeo,...
    segment_same_length,GDPPrior,...
    fraction_change);
J = 2;
pGeo = 1/30;
[ X2,K,rhosTrue,zTrue ] = generateCopulaDataLowDimCP( n,J,pGeo,...
    segment_same_length,GDPPrior,...
    fraction_change);

[n,T] = size(X1);

X = zeros(n,T);
w = 0.5;
for t=1:T
    if rand(1)<w
        X(:,t)= X1(:,t);
    else
        X(:,t) = X2(:,t);
    end
end

tau0 = cumsum(K);
[~,T] = size(zTrue);
n_quadrature_points = 1e2;
z = zTrue;
pGeo0 = 0.5;

% save data to struct
params.X = X;
params.z = randn(1,T);
params.GDPPrior = GDPPrior;
params.pGeo = pGeo0;
params.n_quadrature_points = n_quadrature_points;
params.constrain_rho = false;


%%
experiment_dir = 'experiments/mixture_CP';
if ~exist(experiment_dir, 'dir'),mkdir(experiment_dir),end

N = 10;
M = 10000;


a_prior = 1;
b_prior = 1;
tau = [0];

params.a_prior = a_prior;
params.b_prior = b_prior;

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

save(sprintf('%s/mixture_CP_n%i_N%i',experiment_dir,n,N))















