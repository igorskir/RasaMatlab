function g = imnorm(f, param) 
%IMNORM make nolization of image f 
%param might be 'norm1' or 'norm255' only 
f = double(f); 
f = f - min(f(:)); 
f = f ./ max(f(:)); 
if strcmp(param, 'norm1') 
g = f; 
elseif strcmp(param, 'norm255') 
g = 255*f; 
else 
error('Bad input param argument'); 
end