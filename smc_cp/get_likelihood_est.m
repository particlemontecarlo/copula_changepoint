function [like_est] = get_likelihood_est(log_W_all)
%GET_LIKELIHOOD_EST Summary of this function goes here
%   Detailed explanation goes here

[T,~] = size(log_W_all);

log_Z = zeros(1,T);
for t=1:T
    log_W = log_W_all(t,:);
    log_W_pos = log_W(log_W>-Inf);
    log_Z(t) =  log(sum(exp(log_W_pos-max(log_W_pos)))) + max(log_W_pos);
end

like_est = sum(log_Z);
end

