function d = getBhattacharyyaDistance(X1,X2)
% BHATTACHARYYA  Bhattacharyya distance between two Gaussian classes
%
% d = bhattacharyya(X1,X2) returns the Bhattacharyya distance between X1 and X2.
%
% Inputs: X1 and X2 are n x m matrices represent two sets which have n
%         samples and m variables.
%
% Output: d is the Bhattacharyya distance between these two sets of data.

% Reference:
% Kailath, T., The Divergence and Bhattacharyya Distance Measures in Signal
% Selection, IEEE Trasnactions on Communication Technology, Vol. 15, No. 1,
% pp. 52-60, 1967

%Check inputs and output
error(nargchk(2,2,nargin));
error(nargoutchk(0,1,nargout));

[n,m]=size(X1);
% check dimension 
% assert(isequal(size(X2),[n m]),'Dimension of X1 and X2 mismatch.');
assert(size(X2,2)==m,'Dimension of X1 and X2 mismatch.');

mu1=mean(X1);
C1=cov(X1);
mu2=mean(X2);
C2=cov(X2);
C=(C1+C2)/2;
dmu=(mu1-mu2)/chol(C);
try
    d=0.125*dmu*dmu'+0.5*log(det(C/chol(C1*C2)));
%     d = 0.125*(mu1 - mu2)'*(C^-1)*(mu1 - mu2)+0.5*log(det(C)/sqrtm(det(C1)*det(C2)));
catch
    d=0.125*dmu*dmu'+0.5*log(abs(det(C/sqrtm(C1*C2))));
    warning('MATLAB:divideByZero','Data are almost linear dependent. The results may not be accurate.');
end
% d=0.125*dmu*dmu'+0.25*log(det((C1+C2)/2)^2/(det(C1)*det(C2)));
