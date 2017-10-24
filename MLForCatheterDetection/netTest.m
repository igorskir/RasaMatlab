tic;
outputNet = ones(2, numel(simTargets), 'double');
outputNet(1, :) = simTargets; 
for i = 1:numel(simTargets)
   temp = simInputs(1:end,i);
   outputNet(2,i) = round(net(temp)); 
end
toc;
vars.netTest = {'temp', 'i'};
clear(vars.netTest{:});
outputNet = outputNet';

