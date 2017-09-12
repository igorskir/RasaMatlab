function av = myAverage(A)
% AVERAGE Computes the average value of an array. 
% AV = AVERAGE (A) computes the average value of input 
% array, A, which must be a 1-D or 2-D array. 
% Check the validity of the input. (Keep in mind that 
% a 1-D array is a special case of a 2-D array.) 
if (ndims(A) > 2)
        error('The dimension of the input cannot exceed 2/')
end

av = sum(A(:))/length(A(:));