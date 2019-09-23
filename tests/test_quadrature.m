clear
rng(1)
addpath(genpath('util'))

n = 1;
J = 1;
pGeo = 1;
segment_same_length = true;
GDPPrior = true;

[ X,K,rhos,z ] = generateCopulaData( n,J,pGeo,segment_same_length,GDPPrior );

n_quadrature_points = 1e6;
x = linspace(-1,1);
I = @(x) 5*x.^3-1*x.^2+x+10;
plot(x,I(x))
result = integral(I,-1,1)


thresh = 1e-10;
a = -1+thresh;
b = 1-thresh;
h = (b-a)/(n_quadrature_points);
SimpsonFactor =  [1 repmat([4 2],1,n_quadrature_points/2 - 1) 4 1]';
rhoGrid = linspace(a,b,n_quadrature_points+1)';
pnew = h/3 * sum(I(rhoGrid).*SimpsonFactor)

% use MATLAB's quadrature
like = @(rho) normpdf(X(1),rho*z(1),sqrt(1-(rho.^2)));
prior = @(rho)  GDPdensity(rho./(1-rho.^2).^0.5,3,1/3) .* abs((1 - rho.^2).^(-3/2));
I = @(rho) like(rho).*prior(rho);
plot(x,I(x))
log_int_quad = log(integral(I,-1,1))


% use simspon's rule
pnew = log(h/3 * sum(I(rhoGrid).*SimpsonFactor))





