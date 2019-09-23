function [SS_all,log_W_all,log_W_bar_all] = forwardFilteringFearnhead(params,T)
%%% performs forward filtering recursively. Should return the weights of
%%% the filtering distribution and the support 
SS_0 = [0];
log_W0 = log([1]);
log_Wbar0 = log_gn(1,0,params);

SS_all = zeros(T,T);
SS_all(1,1) = SS_0;
log_W_all = -Inf*ones(T,T);
log_W_all(1,1) = log_W0;
log_W_bar_all = -Inf*ones(T,T);
log_W_bar_all(1,1) = log_Wbar0;

SS_nm1 = SS_0;
log_W_nm1 = log_W0;
for n=2:T
    [SS_updated,log_Wbar_updated] = filteringRecursion(n,SS_nm1,log_W_nm1,params);
    log_W_nm1 = log_Wbar_updated - max(log_Wbar_updated)- ...
            log(sum(exp(log_Wbar_updated-max(log_Wbar_updated))));
    SS_nm1 = SS_updated;
    SS_all(n,1:n) = SS_updated;
    log_W_all(n,1:n) = log_W_nm1;
    log_W_bar_all(n,1:n) = log_Wbar_updated;
end

end


function [SS_updated,log_Wbar_updated] = filteringRecursion(n,SS_nm1,log_W_nm1,params)
%%%% Forward filtering for changepoints

% % % description of variables
% wbar  : 1 x n    weights of random support afte resampling
% SS    : 1 x n    random support vector

[~,card_SS] = size(SS_nm1);
assert(abs(sum(exp(log_W_nm1))-1)<1e-8)

% get value of last Wbar
summand_arr = zeros(1,card_SS);
for i=1:card_SS
    xnm1 = SS_nm1(i);
    log_W_val = log_W_nm1(i);
    summand_arr(i) = log_fn(n-1,xnm1,params) + log_W_val;
end

log_Wbar_n_nm1 = log_gn(n,n-1,params)+ ...
        log(sum(exp(summand_arr-max(summand_arr))))+max(summand_arr);


% get values of other Wbars
log_Wbar_arr = zeros(1,card_SS);
for i=1:card_SS
    xnm1 = SS_nm1(i);
    log_W_val = log_W_nm1(i);
    log_Wbar_arr(i) = log_gn(n,xnm1,params) +...
        log_fn(xnm1,xnm1,params) + log_W_val;      
end

% update random support
cardSS_n = n;
SS_updated = zeros(1,cardSS_n);
SS_updated(1:(cardSS_n-1)) = SS_nm1;
SS_updated(cardSS_n) = n-1;

% update weights of random support
log_Wbar_updated = zeros(1,cardSS_n);
log_Wbar_updated(1:(cardSS_n-1)) = log_Wbar_arr;
log_Wbar_updated(cardSS_n) = log_Wbar_n_nm1;

end







