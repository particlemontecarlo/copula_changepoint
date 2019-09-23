function [n_changepoints] = n_changepoints_from_list(CPall)
%N_CHANGEPOINTS_FROM_LIST Counts the number of changepoints in each sample

f = find(diff([0,CPall>0,0]==1));
p = f(1:2:end-1);  % Start indices
n_changepoints = f(2:2:end)-p;  % Consecutive ones? counts

end

