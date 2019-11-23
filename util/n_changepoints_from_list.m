function [n_changepoints] = n_changepoints_from_list(CPall)
%N_CHANGEPOINTS_FROM_LIST Counts the number of changepoints in each sample

if ~any(diff(CPall)==0) && CPall(end)~=0
    f = find(diff([0,CPall>0,0]==1));
    p = f(1:2:end-1);  % Start indices
    n_changepoints = f(2:2:end)-p;  % Consecutive ones? counts
else
    n_changepoints = get_num_cp(CPall);
end
end

function [segment_length] = get_num_cp(tau_collect)
  
    segment_length = zeros(1,sum(tau_collect==0));
    
    zero_indices = find(tau_collect==0);
    for i=1:length(zero_indices)
        if zero_indices(i)~=zero_indices(end)
            start_indx = zero_indices(i);
            end_indx = zero_indices(i+1);
            segment_length(i) = end_indx-start_indx-1;
        else
            segment_length(i) = length(tau_collect)-zero_indices(end);
        end
        
    end
end





