numFeats = numel(featureNames);
% numFeats = 2;
score = zeros(1, numFeats-1);
j = 1;
for i = 2:numFeats
    CreateProbabilityPlot(dataCatheter(:,i),dataTissue(:,i), featureNames{1,i});
    score(1,j)= GetDensityDifference(dataCatheter(:,i),dataTissue(:,i), 'ShowPlot');
    j = j+1;
    pause(5)
    delete(findall(0,'Type','figure'))
    bdclose('all')
end
