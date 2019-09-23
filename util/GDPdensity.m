function [ p,varargout  ] = GDPdensity( lambda,alpha,xi )
%GDPDENSITY Evaluates density of generalised double Pareto prior, as
%recommended in Murray (2011)


p = 1/(2*xi) * ( 1+abs(lambda)/(alpha*xi) ).^-(alpha+1);

if nargout==2
    varargout{1} = log(p);
end

end

