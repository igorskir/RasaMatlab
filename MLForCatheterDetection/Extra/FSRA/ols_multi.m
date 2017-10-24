function [r,beta]=ols_multi(y,x,xx,exc)

% OLS_MULTI - Performs a multivariate OLS regression.
%
% Note: 'x' shall not contain a constant!
%
% Author: Marco Buchmann (mb-econ.net), ECB / Frankfurt University
% Date: March 2009
%
% Required input:
% -- Tx1 vector 'y' that holds data for the explained variable
% -- TxK matrix 'x' that holds data for K explanatory variable(s)
% -- SxK matrix 'xx' at which to estimate the regression function
%
% Optional input:
% -- scalar 'exc': if exc=1, then constant is excluded from the model 
%
% Output:
% -- Sx1 vector 'r', the function fit at point/vector/matrix 'xx'
% -- 1x(K+1) vector 'beta' that holds the constant and derivative estimates

% Check input
if nargin<3
    error('Not enough input arguments')
end
if size(y,1)~=size(x,1)
    error('Dimension mismatch b/w dependent and independent variables');
end
if size(x,2)~=size(xx,2)
    error('Dimension mismatch')
end
if exist('exc')==1 %#ok<EXIST>
    ex=exc;
else
    ex=0;
end

% Estimation + generate fit
if ex==1 % then exclude constant
    beta=(x'*x)\x'*y;
    r=xx*beta;
else     % include constant  
    xz=[ones(size(x,1),1) x];
    beta=(xz'*xz)\xz'*y;
    r=[ones(size(xx,1),1) xx]*beta;
end