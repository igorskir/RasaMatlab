[M,N] = size(A);

ColSum = zeros(M,1);
RowSum = zeros(N,1);
A = double(A);
tic();
for i = 1:M
    for j = 1:N
        ColSum(j) = ColSum(j) + A(i,j);
        RowSum(i) = RowSum(i) + A(i,j);
    end
end
toc();

%% New Alorithms
T = A;
AM = T(:); 
for i = 1:M*N
    
    if (AM(i) < 68)
        AM(i) = 0;
        continue;
    end
    AM(i) = 255;
    
end

T(:) = AM(:); 