function av = average(A)
%AVERAGE Computes the average value of an array. 
if (ndims(A) > 2)
    error('cannot exceed 2/');
end
av = sum(A(:))/length(A(:));