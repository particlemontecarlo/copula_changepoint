clear

tau_collect = [0,0,2,4,0,2,0,2,0];
T = 6;
[bin_arr] = binary_changepoints_from_list(tau_collect,T)