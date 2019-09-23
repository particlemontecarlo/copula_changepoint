% this captures the predictive likelihood conditional on those in the
% current segment
function [log_g_val] = log_gn(n,xnm1,params)

Yall = params.Y;
yn = Yall(n);
Y_conditioning = Yall((xnm1+1):(n-1));
[~,log_py] = pred_like_est(yn,Y_conditioning,params);

log_g_val = log_py;

end

function [py,log_py] = pred_like_est(yn,Y,params)
%PRED_LIKE_EST predictive likelihood for normal model

sigma02 = params.sigma02;
sigma2 = params.sigma2;

[~,n] = size(Y);

if n>0
    sigma02_dash = (1/sigma02 + n/sigma2)^-1;
    mu0_dash = sigma02_dash*( sum(Y) / sigma2);
else
    sigma02_dash = (1/sigma02 + n/sigma2)^-1;
    mu0_dash = sigma02_dash*( sum(Y) / sigma2);
end

s2_py = sigma02_dash+sigma2;
log_py = -0.5*(log(2*pi*s2_py) + (yn-mu0_dash)^2/s2_py);
py = exp(log_py);
end

