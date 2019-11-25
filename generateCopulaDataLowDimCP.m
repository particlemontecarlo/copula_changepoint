function [ X,K,rhos,zTrue ] = generateCopulaDataLowDimCP( n,J,pGeo,...
    segment_same_length,GDPPrior,...
    fraction_change)
%GENERATEDATA want a method that generates copula data, though where the
%changes in the factor loadings only occurs in a small number of components
%of the total number of latent factor loadings


% get the number of components to change
change_components = ceil(fraction_change*n);

K = 0;
if ~segment_same_length
    for j=1:J
        K = [K geornd(pGeo)+1];
    end
else
    for j=1:J
        K = [K round(1/pGeo)];
    end
end

% K = [0 L*ones(1,J)];
% p = J/sum(K);
T = sum(K);


%Sample latent factors
z = (randn(1,T));
%z = 1*(randn(1,T)-0.5);

% Sample conditionally on segment
X = zeros(n,T);
rhos = zeros(n,J);
for j = 1:J
    
    if GDPPrior
        lambda = sampleGDP( 3,1/3,n );
    else
        lambda = randn(n,1);
        warning('Generating from normal prior')
    end
    
    rho_new =  lambda./(1+lambda.^2).^0.5;
    
    if j==1
        rho = rho_new;
    else
        rho(1:change_components) = rho_new(1:change_components);
    end
    rhos(:,j) = rho;

    
    tstart = sum(K(1:j))+1;
    tend = sum(K(1:(j+1)));
    for t = tstart:tend
        for i = 1:n
             X(i,t) = normrnd(rho(i)*z(t),(1-rho(i)^2)^0.5);
        end
    end
end
zTrue = z;
end

