function misclass = svmFun(xTrain, yTrain, xTest, yTest) 

x = input;
y = output;
crossValidation = cvpartition(y, 'Holdout', 0.15);
xTrain = x(crossValidation.training,:);
yTrain = y(crossValidation.training);
xTest = x(crossValidation.test,:);
yTest = y(crossValidation.test);
%%
tic
% SVMModel = fitcsvm(xTrain, yTrain);
SVMModel = fitcsvm(xTrain, yTrain,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');
CVSVMModel = crossval(SVMModel);
classLoss = kfoldLoss(CVSVMModel)
% KNNModel = fitcknn(xTrain, yTrain, 'OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
% KNNModel = fitcknn(xTrain, yTrain, 'Distance', 'seuclidean');
toc

%%
yPred = predict(SVMModel, xTest);
misclass = sum(yTest ~= yPred);
end