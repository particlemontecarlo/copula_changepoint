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
J = 2;
pGeo = 0.02;
segment_same_length = true;
GDPPrior = false;
n_quadrature_points = 1e2;
constrain_rho = false;


[ X,K,rhosTrue,zTrue ] = generateCopulaData( n,J,pGeo,segment_same_length,GDPPrior );
corr_mat_true = rhosTrue(:,1)*rhosTrue(:,1)';



pGeo0 = 0.5;

[~,T] = size(zTrue);
tau0 = [0,1:10:T];

% save data to struct
params.X = X;
params.GDPPrior = GDPPrior;
params.pGeo = pGeo0;
params.n_quadrature_points = n_quadrature_points;



%%
n_vals = [10];
experiment_dir = 'experiments/compare_switching';
if ~exist(experiment_dir, 'dir'),mkdir(experiment_dir),end

N = 10;
M = 1000;
n_M = 10;
constrain_rho_arr = [false,true];


for ss=1:length(constrain_rho_arr)
    constrain_rho = constrain_rho_arr(ss);
    
    for n_i=1:n_M
        
        a_prior = 1;
        b_prior = 1;
        tau = tau0;
        params.z = randn(1,T);
        params.a_prior = a_prior;
        params.b_prior = b_prior;
        params.constrain_rho = constrain_rho;

        tau_collect = [tau];
        like_ests = zeros(1,M);
        pGeo_collect = zeros(1,M);
        corr_mat_collect =  zeros(n,n,M);
        z_collect = zeros(M,T);
        joint_loglik_collect = zeros(1,M);
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
            joint_loglik_collect(m) = joint_loglik;
            rho_record_collect(m,:) = rho_record;
        end
        time_taken = toc;

        save(sprintf('%s/compare_switching_n%i_constrain%s',...
            experiment_dir,n_i,string(constrain_rho)))
    end


end











