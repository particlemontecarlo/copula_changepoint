function [ O_nm1 ] = stratifiedResampling( I_nm1,W_for_I_nm1,N,L_nm1 )
%STRATIFIEDRESAMPLING Resamples particles with weights below a threshold
%   As the space is discrete the resampling corresponds only to 
%   reweighting of the particles
%   I_nm1 is the location of 'unsafe' particles - i.e. with low weight
%   L_nm1 is the number of 'safe' particles - i.e. with high weight
W_nm1 = W_for_I_nm1;

% check the location of the unsafe particles corresponds to the 
assert(all(size(I_nm1)==size(W_nm1)));

% check that the total number of particlces is correct
[~,n_not_safe]= size(I_nm1);
assert((n_not_safe+L_nm1)<=(N+1)) 


% normalise weights
Wbar_nm1 = W_nm1/sum(W_nm1);

% construct CDF
Q_nm1 = [0,cumsum(Wbar_nm1)];

% generate u variables
U1 = rand(1)/(N-L_nm1);
U_arr = U1 +  (0:(N-L_nm1-1))/(N-L_nm1);

% set off spring to 0 or 1
O_nm1 = histcounts(U_arr,Q_nm1);


assert(all(O_nm1==0|O_nm1==1))

end










