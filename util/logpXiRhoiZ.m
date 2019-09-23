function [ logp ] = logpXiRhoiZ( X,z,rho,GDPPrior )

prior_func = @(Rho) priorFunc(X,z,Rho,GDPPrior);
logp = arrayfun(prior_func,rho);
end

function logp = priorFunc(X,z,rho,GDPPrior)

%Jacobian  = rho.^2./(1 - rho.^2).^(3/2) + (1 - rho.^2).^(-1/2) =
%(1-rho^2)^(-3/2)
Var = (1-rho^2);
logVar = log(Var);


if(length(rho)~=1)
    error('Non-scalar rho in density evaluation');
end

% this ensures that if rho should be positive then logp is not -Inf
logic_condition = (rho>=-1) && (rho<=1);

% evaluate 
if(logic_condition )
    logJacobian = -1.5*logVar;
    if GDPPrior
        logPrior = log(GDPdensity( rho/((1-rho^2)^0.5),3,1/3 ));
    else
        logPrior = log(normpdf(rho/(1-rho^2)^0.5));
    end
    logp =  lognormpdfXRho(X,rho*z,Var,logVar) + logPrior + logJacobian;
else
    logp = -Inf;
end


end

function logp = lognormpdfXRho(X,mu,Var,logVar)

T = length(X);


log1 =  -(T/2)*log(2*pi);
log2 = - (T/2)*logVar;
log3 = - (1/(2*Var))*sum((X-mu).^2);
logp =  log1+log2+log3;

end