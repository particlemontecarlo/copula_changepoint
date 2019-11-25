function [bin_arr] = binary_changepoints_from_list(tau_collect,T)
%binary_changepoints_from_list Converts [0] separated list to array
%indicating where changepoints occurred

n_iter = sum(tau_collect==0);

bin_arr = zeros(n_iter,T);
segment_length = zeros(1,n_iter);

zero_indices = find(tau_collect==0);
for i=1:length(zero_indices)
    
    zero_line = zeros(1,T);
    
    if zero_indices(i)~=zero_indices(end)
        start_indx = zero_indices(i);
        end_indx = zero_indices(i+1);
        
    else
        start_indx = zero_indices(end);
        end_indx = length(tau_collect)+1;
    end
    
    tau_line = tau_collect((start_indx+1):(end_indx-1));
    zero_line(tau_line)=true;
    bin_arr(i,:) = zero_line;
end


end





