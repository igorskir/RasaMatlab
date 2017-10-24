%% Script: Illustrate the use of the FSR Algorithm

clear all; clc;

% Settings
N=100;                   % number of observations
beta=[0.6 0.8 0.3 0.1];  % true regression coefficients

% DGP
x=randn(N,size(beta,2)); % covariates normally distributed
y=x*beta';               

% Call FSRA: 'out' will be a row vector with numbers refering to the order
% of variables as selected by the FSR algorithm
out=fsra(y,x);
