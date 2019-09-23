function [SS_all,log_W_all,log_W_bar_all] = forwardFilteringSMC(N,params)
%%% performs Fearnhead (2007) SMC to obtain filtering approximations for
%%% changepoints
%

[~,T] = size(params.Y);
SS_0 = [0];
log_W0 = log([1]);
log_Wbar0 = log_gn(1,0,params);

% intialise random support and starting weights
SS_all = zeros(T,min(N+1,T));
SS_all(1,1) = SS_0;
log_W_all = -Inf*ones(T,min(N+1,T));
log_W_all(1,1) = log_W0;
log_W_bar_all = -Inf*ones(T,T);
log_W_bar_all(1,1) = log_Wbar0;

% iterate through time, recording random support and weights
SS_nm1 = SS_0;
log_W_nm1 = log_W0;
for n=2:T
    [SS_updated,log_Wbar_updated] = forwardFilteringSMCRecursion(n,N,SS_nm1,log_W_nm1,params);
    log_W_nm1 = log_Wbar_updated - max(log_Wbar_updated)- ...
            log(sum(exp(log_Wbar_updated-max(log_Wbar_updated))));
    SS_nm1 = SS_updated;
    length_SS = min(N+1,n);
    SS_all(n,1:length_SS) = SS_updated;
    log_W_all(n,1:length_SS) = log_W_nm1;
    log_W_bar_all(n,1:length_SS) = log_Wbar_updated;
end

end
