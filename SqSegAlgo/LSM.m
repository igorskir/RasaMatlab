function [ k] = LSM( y )
% Least Square Method
% returns k coeff (y=kx+b) and value of functional
x=1:5;
n=5;
sumX=15;
sumXSq=55;
k=(n*sum(x.*y)-sumX*sum(y))/(n*sumXSq-225);
%k=(n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.*x)-(sum(x))^2);
%b=(sum(y)-k*sum(x))/n;
%Y=k.*x+b;
%F=sum((k.*x+b-y).^2);
end

