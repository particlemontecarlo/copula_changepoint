%% script to generate and plot synthetic data
clear
rng(1)
addpath(genpath('util'))

n = 10;
J = 5;
pGeo = 0.02;
segment_same_length = true;
GDPPrior = true;

[ X,K,rhos,zTrue ] = generateCopulaData( n,J,pGeo,segment_same_length,GDPPrior );

plot(X(1,:))


