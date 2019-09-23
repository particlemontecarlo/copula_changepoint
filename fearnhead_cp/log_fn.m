% we specialise for the geometric distribution
function [log_f_val] = log_fn(xn,xnm1,theta)

if xn==xnm1
    f_val = 1-theta.pGeo;
else
    f_val = theta.pGeo;
end
log_f_val = log(f_val);


end