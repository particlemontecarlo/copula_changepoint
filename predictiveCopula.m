function [ log_g_val ] = predictiveCopula( n,xnm1, params )
%CONDITIONALSMC Performs conditional SMC sampling of changepoints,
%conditioned on an existing changepoint configuration
%   Method follows that of Whiteley (2010) and Fearnhead (2007). Method
%   assumes there exists a changepoint at t=0, for a series indexed between
%   1:t

% unpack arguments
X = params.X;
z = params.z;
GDPPrior = params.GDPPrior;
constrain_rho = params.constrain_rho;
n_quadrature_points = params.n_quadrature_points;

[~,T] = size(X);

conditioning_indices = (xnm1+1):(n-1);

if ~isempty(conditioning_indices)
    for t=conditioning_indices
        if t==conditioning_indices(1)
            % use first index to initialise static data
            [ ~,dynamic_data,static_data ] = predictive_start( ...
                            X(:,t),z(t),GDPPrior,n_quadrature_points,constrain_rho );
        else
            % estimate recursively p(y_t|y_s:tm1)
            pre_calc_t = pre_calculation(static_data,X(:,t),z(t),constrain_rho);
            [ ~,dynamic_data] = predictive(dynamic_data,static_data,pre_calc_t );
                                
        end
    end
    % get the predictive using the final index
    assert((t+1)==n)
    pre_calc_t = pre_calculation(static_data,X(:,n),z(n),constrain_rho);
    [log_g_val,~] = predictive(dynamic_data,static_data,pre_calc_t );
else
    [ log_g_val,~,~ ] = predictive_start( ...
        X(:,n),z(n),GDPPrior,n_quadrature_points,constrain_rho );
end


end


% this can be used as the additional density evaluations need only be
% evaluated once for all the different particles
function [pre_calc] = pre_calculation(static_data,Xt,zt,constrain_rho)

rhoGrid = static_data.rhoGrid;
rhoVar = static_data.rhoVar;
logf1 = logpXcondRhoZ(rhoGrid,Xt',zt,rhoVar);

% apply constraint to integral
if constrain_rho
    [nx,ny] = size(logf1);
    zero_prior_mask = repelem(false,nx,ny);
    zero_prior_mask(1:floor(nx/2),1) = 1;
    logf1(zero_prior_mask) = -Inf;
end
pre_calc = logf1;

end


%Recursively evaluates the posterior predictive
function [logp,dynamic_data] = predictive(dynamic_data,static_data,pre_calc)


logfProdOld = dynamic_data.logfProd;

h = static_data.h;
logrhoInt = static_data.logrhoInt;
logdenom = dynamic_data.logdenom;


%logf1 = logpXcondRhoZ(rhoGrid,Xt',zt,rhoVar);
logf1 = pre_calc;
logfProdNew = logf1 + logfProdOld;
a1 = logfProdNew + logrhoInt - logdenom;
max_a1 = max(a1,[],1);
diff_mat = bsxfun(@minus,a1,max_a1);
logpAll = log(h/3) + log(sum(exp(diff_mat ))) + max_a1 ;

dynamic_data.logfProd = logfProdNew;
dynamic_data.logdenom = logpAll + logdenom;

if any(isinf(logpAll))
    warning('Underflow in prediction')
    logpAll = logpAll(~isinf(logpAll));
end


logp = sum(logpAll);

end


function [ logp,dynamic_data,static_data ] = predictive_start( ...
    Xt,zt,GDPPrior,n_quadrature_points,constrain_rho )
% PREDICTIVESTART sets up the quadrature rules and the static data settings

thresh = 1e-10;
a = -1+thresh;
b = 1-thresh;


h = (b-a)/(n_quadrature_points);
SimpsonFactor =  [1 repmat([4 2],1,n_quadrature_points/2 - 1) 4 1]';
rhoGrid = linspace(a,b,n_quadrature_points+1)';

if GDPPrior
    rhoPrior = rhoPriorGDP(rhoGrid);
else
    rhoPrior = rhoPriorNormal(rhoGrid);
end
logrhoInt = log(rhoPrior)+log(SimpsonFactor);
rhoVar = (1-rhoGrid.^2);


logf1 = logpXcondRhoZ(rhoGrid,Xt',zt,rhoVar);
% apply constraint to integral
if constrain_rho
    [nx,ny] = size(logf1);
    zero_prior_mask = repelem(false,nx,ny);
    zero_prior_mask(1:floor(nx/2),1) = true;
    logf1(zero_prior_mask) = -Inf;
end

pnew = h/3 * sum(exp(logf1+logrhoInt),1);
logfProd = logf1;
logdenom = log(pnew);

logp = sum(logdenom);

dynamic_data.logfProd = logfProd;
dynamic_data.logdenom = logdenom;

static_data.rhoGrid = rhoGrid;
static_data.logrhoInt = logrhoInt;
static_data.rhoVar = rhoVar;
static_data.h = h;

end



function logp = logpXcondRhoZ(Rho,Xt,zt,rhoVar)

mu = Rho*zt;
logp = - 0.5*log(2*pi*rhoVar) - 0.5*((Xt-mu).^2)./rhoVar ;

end



%Generalised double Pareto density
function rhoPrior = rhoPriorGDP(Rho)

rhoPrior = GDPdensity(Rho./(1-Rho.^2).^0.5,3,1/3) .* abs((1 - Rho.^2).^(-3/2));
end

%Normal prior density
function rhoPrior = rhoPriorNormal(Rho)

rhoPrior = normpdf(Rho./(1-Rho.^2).^0.5) .* abs((1 - Rho.^2).^(-3/2));
end


%Generalised Pareto density for first component
function rhoPrior = rhoPriorGP(Rho)

rhoPrior = GPdensity(Rho./(1-Rho.^2).^0.5,3,1/3) .* abs((1 - Rho.^2).^(-3/2));
end



