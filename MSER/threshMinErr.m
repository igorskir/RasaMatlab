%function [t,s]=threshMinErr(h,x=1:length(h))
%
% Kittler & Illingworth - Minimum error threshold.
%
%Input:
% h - histogram of data
% x - histogram x-axis (default is 1:length(h))
%
%Output:
% t - suggested threshold
% s - the error function that is minimized
%
%
% Calculates the threshold value needed to threshold the data h
% according to the minimum error thresholding metods as described in 
% J. Kittler & J. Illingworth: "Minimum Error Thresholding"
%  Pattern Recognition, Vol 19, nr 1. 1986, pp. 41-47.
%
% Algorithm is the direct (non-iterative) version from
% Glasbey, C. A. 1993. "An analysis of histogram-based thresholding algorithms."
%  Graphical Models and Image Processing, 55(6): 532-7. 
%
%
% Joakim Lindblad 2002-08-26


function [t,s]=threshMinErr(h,x)

h=h(:);
n=length(h);

if nargin<2
	x=1:length(h);
end
x=x(:);

A=cumsum(h);
B=cumsum(h.*x);
C=cumsum(h.*x.^2);

warning off %zero divides

p=A/A(n);
q=(A(n)-A)/A(n);

u=B./A;
v=(B(n)-B)./(A(n)-A);

s2=C./A-u.^2;
t2=(C(n)-C)./(A(n)-A)-v.^2;

s=( p.*log(s2./p) + q.*log(t2./q) );

warning on

f=find(isfinite(s));
[dummy,idx]=min(s(f));
t=x(f(idx));