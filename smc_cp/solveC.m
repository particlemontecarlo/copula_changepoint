% Finds C, through solving sum(min(1,w_i C )) = N
% Requires sorting W and finding the solution through the solving N =
% min(1,C*w_i) as the sum of piecewise linear equations
function [C,logC] = solveC(logW)

% check weights are normalised
if abs(sum(exp(logW))-1)>1e-8
    error('Weights must sum to 1')
end
 % make sure all weights are positive
if any(logW==-Inf)
    error('Particle weights must be positive')
end

% get the number of particles
[~,Np1] = size(logW);
assert(Np1>1)

% number of particles
N = Np1-1;


W = exp(logW);

% want sorted largest to smallest
W_sorted = sort(W,'descend');

x_vals = 1./W_sorted;
g_vals = cumsum(W_sorted,'reverse');
y_vals = zeros(1,Np1);
y_vals(1) = g_vals(1)*x_vals(1);

ii = 2;
xstar = Inf;

% test first one
if y_vals(1)>N
    xstar = N/g_vals(1);
else
    if N>1
        % test rest
        while ii<=Np1 && xstar==Inf
            
            y_vals(ii) = y_vals(ii-1)+(x_vals(ii)-x_vals(ii-1))*g_vals(ii);
            
            if y_vals(ii)>N
                xstar = (N-y_vals(ii-1))/g_vals(ii) + x_vals(ii-1);
            end
            ii = ii+1;
        end
    else
        error('Resampling one particle')
    end
end


C = xstar;
logC = log(C);
Ntest = sum(min(1,W*C));

abs_err = abs(Ntest-N);
if abs_err>1e-8
    warning('abs err in resampling exceeding 1e-8')
end

end



% 
% % Finds C, through solving sum(min(1,w_i C )) = N
% % Requires sorting W and finding the solution through the solving N =
% % min(1,C*w_i) as the sum of piecewise linear equations
% function [C,logC] = solveC_dev1(logW)
% 
% % check weights are normalised
% assert(sum(exp(logW))==1)
% 
% % get the number of particles
% [~,Np1] = size(logW);
% assert(Np1>1)
% 
% % number of particles
% N = Np1-1;
% 
% 
% W = exp(logW);
% 
% % want sorted largest to smallest
% W_sorted = sort(W,'descend');
% 
% x_vals = 1./W_sorted;
% g_vals = cumsum(W_sorted,'reverse');
% y_vals = zeros(1,Np1);
% 
% ii = N;
% xstar = Inf;
% y_vals(Np1) = Np1; 
% while xstar==Inf
%     y_vals(ii) = y_vals(ii+1) - g_vals(ii+1)*(x_vals(ii+1)-x_vals(ii));
%     if y_vals(ii)<N && xstar==Inf
%         xstar = (N-y_vals(ii))/g_vals(ii+1) + x_vals(ii);
%     end
%     ii = ii-1;
% end
% 
% C = xstar;
% logC = log(C);
% 
% 
% end
