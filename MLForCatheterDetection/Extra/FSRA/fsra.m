function out=fsra(y,x)

% Forward Stepwise Regression Algorithm (FSRA).
%
% Author: Marco Buchmann (mb-econ.net), ECB / Frankfurt University
% Date: April 2010
%
% Required input:
% -- a Tx1 vector 'y': the dependent variable 
% -- a TxK matrix 'x': the set of explanatory variables (no constant to be included)
%
% Output:
% -- Kx1 vector holding integers that refer to the variables in 'x'

% Check input
if nargin<2
   error('Not enough input arguments') 
end
if size(y,1)~=size(x,1)
    error('Number of observations for dependent and independent variables does not match')
end
tmpv=var(x); tmpt=(tmpv==0);
if sum(tmpt)>0
    error('Please exclude the constant from X')
end
if size(x,2)<2
    error('At least two regressors to be provided')
end

% Set parameters and pre-allocate space for output
T=size(y,1); K=size(x,2); out=NaN(1,K);

% Normalize 'y' and 'x'
y=y-mean(y);
for k=1:K
   x(:,k)=(x(:,k)-mean(x(:,k)))./std(x(:,k)); 
end
xx=x; % duplicate 'x'

% FSR Algorithm
tar=y;
for i=1:K-1 % outer loop: over covariates
    KK=size(xx,2);                                % no. of columns of current 'xx' matrix
    beta=NaN(KK,1);                               % space for coefficient estimates
    resid=NaN(T,KK);                              % space for residuals
    if i>1 
        [r1,r2]=ols_multi(y,x(:,out(1,1:i-1)),x(:,out(1,1:i-1)),1); % regression of 'y' on latest rel. set of covariates  
        tar=y-r1;                                 % residuals  
    end
    for k=1:KK                                    % loop over covariates
        [t1,t2]=ols_multi(tar,xx(:,k),xx(:,k),1); % projection, excluding constant
        resid(:,k)=tar-t1;                        % residuals
        beta(k,1)=t2;                             % save coef. estimate
    end
    [t1,t2]=max(abs(beta));                       % find max absolute correlation
    tmpd=NaN(T,K);                                % determine correct identifier for covariate
    for u=1:K
       tmpd(:,u)=x(:,u)-xx(:,t2); 
    end
    sumx=sum(tmpd,1);   
    [s1,s2]=max((sumx==0));
    out(1,i)=s2;                                  % save intermediate result
    xx(:,t2)=[];                                  % exclude latest rel. covariate
    tar=resid(:,t2);                              % new target
end
out(1,end)=sum(linspace(1,K,K))-nansum(out);