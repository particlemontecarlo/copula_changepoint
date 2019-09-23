function [SS_updated,log_Wbar_updated] = forwardFilteringSMCRecursion(n,N,SS_nm1,log_W_nm1,params)
%%%% Forward filtering for changepoints

% % % description of variables
% wbar  : 1 x min(N+1,n)    weights of random support afte resampling
% SS    : 1 x min(N+1,n)    random support vector
% S     : 1 x min(N+1,n)    indicator array
% In    : 1 x min(N+1,n)    array used in resampling


% we have that SS captures where a possible changepoint is with SSn at most 
% roughly of length n and that Sn is an indicator array of size n, 
% suggesting whether the particle is alive or dead
% for each particle indexed by elements in SSn we have a weight
[~,card_SS] = size(SS_nm1);

%%%% need
% 1) deal properly when some input weights can be zero
% 2) assign proper weights to those that resampled - this is dealt with
% with the min in
assert(abs(sum(exp(log_W_nm1))-1)<1e-8)

if n<=(N+1)
    % no need to resample for small values of n
    C = Inf;
    logCnm1 = Inf;
    S_nm1 = ones(1,n-1);
else
    % resample if N is large enough
    % note that this requires only positive weights to make sense
    [C,logCnm1] = solveC(log_W_nm1(log_W_nm1~=-Inf));
    
    % get the unsafe particles
    unsafe_mask = log_W_nm1+logCnm1<=0;
    I_nm1 = SS_nm1(unsafe_mask);
    L_nm1 = card_SS-sum(unsafe_mask);
    W_for_I_nm1 = exp(log_W_nm1(unsafe_mask));
    assert(sum(unsafe_mask)>0)

    [ O_nm1 ] = stratifiedResampling( I_nm1,W_for_I_nm1,N,L_nm1 );

    S_nm1 = zeros(1,N);
    S_nm1(~unsafe_mask) = 1;
    S_nm1(unsafe_mask) = O_nm1;
    assert(all((S_nm1==0)|(S_nm1==1)));
end

% convert S_nm1 to logical array
S_nm1 = S_nm1&1;


% get value of last Wbar
summand_arr = zeros(1,card_SS);
for i=1:card_SS
    xnm1 = SS_nm1(i);
    log_W_val = log_W_nm1(i);
    summand_arr(i) = log_fn(n-1,xnm1,params) + log_W_val;
end

%Wbar_n_nm1 = gn(n-1,n-1,params)sum(summand_arr);
log_Wbar_n_nm1 = log_gn(n,n-1,params)+ ...
        log(sum(exp(summand_arr-max(summand_arr))))+max(summand_arr);

% get values of other Wbars
log_Wbar_arr = zeros(1,card_SS);
for i=1:card_SS
    xnm1 = SS_nm1(i);
    log_W_val = log_W_nm1(i);
    log_Wbar_arr(i) = log_gn(n,xnm1,params) +...
        log_fn(xnm1,xnm1,params) + log_W_val - ...
        min(0,log_W_val + logCnm1);      
end

% update random support
cardSS_n = min(N+1,n);
SS_updated = zeros(1,cardSS_n);
SS_updated(1:(cardSS_n-1)) = SS_nm1(S_nm1);
SS_updated(cardSS_n) = n-1;

% update weights of random support
log_Wbar_updated = zeros(1,cardSS_n);
log_Wbar_updated(1:(cardSS_n-1)) = log_Wbar_arr(S_nm1);
log_Wbar_updated(cardSS_n) = log_Wbar_n_nm1;

end







