%% estimates posterior predictive
% the following script estimates the posterior predictive for observations
% arising in the factor copula model. The methodology is the start of 
clear
rng(1)
addpath(genpath('util'))

n = 5;
J = 1;
pGeo = 0.5;
segment_same_length = true;
GDPPrior = true;

[ X,K,rhos,zTrue ] = generateCopulaData( n,J,pGeo,segment_same_length,GDPPrior );

n_quadrature_points = 1e3;
z = zTrue;


params.X = X;
params.z = z;
params.GDPPrior = GDPPrior;
params.n_quadrature_points = n_quadrature_points;
[ loglik_record ] = predictiveCopula( 2,0,params );
loglik_record




%% the integral that is being performed is
% for i in 1:n
%       int N(X_t^{i}|rho_i*z_t,1-rho_i^2)*pi(rho_i) d(rho_i)
% we test on the first two observations of the generated data, where the
% second is a conditional likelihood of p(X_2|X_1) which we estimate using
% the ratio of p(X_1:2)/p(X_1)
thresh = 1e-10;
a = -1+thresh;
b = 1-thresh;
waypoint_mesh = linspace(a,b,n_quadrature_points+1);

%Generalised double Pareto density
rhoPrior = @(rho)  GDPdensity(rho./(1-rho.^2).^0.5,3,1/3) .* abs((1 - rho.^2).^(-3/2));

int_f_vals = zeros(1,n);
for i=1:n
    X_ti = X(i,1);
    zt = z(1);
    rhoLikelihood = @(rho) normpdf(X_ti,rho*zt,sqrt(1-rho.^2));
    int_f = @(rho) rhoLikelihood(rho).*rhoPrior(rho);
    int_f_vals(i) = integral(int_f,-1,1,'RelTol',0,'AbsTol',1e-10,'waypoints',waypoint_mesh);    
end
log_pX1  = sum(log(int_f_vals))

% perform same calculation for first two observations
int_f_vals = zeros(1,n);
for i=1:n
    rhoLikelihood1 = @(rho) normpdf(X(i,1),rho*z(1),sqrt(1-rho.^2));
    rhoLikelihood2 = @(rho) normpdf(X(i,2),rho*z(2),sqrt(1-rho.^2));
    int_f = @(rho)  rhoLikelihood1(rho).*rhoLikelihood2(rho).*rhoPrior(rho);
    int_f_vals(i) = integral(int_f,-1,1,'RelTol',0,'AbsTol',1e-10,'waypoints',waypoint_mesh);    
end
log_pX2g1  = sum(log(int_f_vals))-log_pX1
















