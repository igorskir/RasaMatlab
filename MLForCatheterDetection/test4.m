alpha = 0.05;
%% Normality
for i = 1:7
    [h, p] = kstest(b(:,i));
end
%% The same distribution
j = 1;
for i = 1:6
    [resWilcoxon(1,j), resWilcoxon(2,j)] = ranksum(b(:,1), b(:,i+1),'alpha',alpha);
    j = j + 1;
end