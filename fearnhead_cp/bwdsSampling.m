function [ tau ] = bwdsSampling(params,SS_all,log_W_all)
%BWDSSAMPLING Performs the perfect sampling of posterior changepoint
%locations given filtering recursions
%   Implements the posterior sampling as given in Fearnhead (2007) and
%   Whiteley (2011)

[T,~] = size(SS_all);


% backward sampling
taui = SS_all(T,find(  mnrnd(1,exp(log_W_all(T,:)))   ));
tau = taui;
while taui>1

    log_Wtau = log_W_all(taui,:);
    Stau = SS_all(taui,:);
    
    logWsupport = log_Wtau(log_Wtau>-Inf);
    Ssupport = Stau(log_Wtau>-Inf);
    nsupport = length(Ssupport);
    
    logpsupport = zeros(1,nsupport);
    for i = 1:nsupport
        logpsupport(i) = logWsupport(i) + log_fn(taui,Ssupport(i),params);
    end
    logpn = logpsupport-max(logpsupport);
    taui = Ssupport(  find(mnrnd(1,exp(logpn)/sum(exp(logpn))))  );
    tau = [taui tau];
end
if taui~=0
    tau = [0 tau];
end


end

