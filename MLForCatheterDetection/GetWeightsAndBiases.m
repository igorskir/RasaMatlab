function [B,IW,LW] = GetWeightsAndBiases(inputs, targets)
    try
        if size(inputs,1) ~= size(targets,1)
            error('Inputs and Targets should have the same lenght')
        end

        if sum(isnan(inputs(:))) > 0
            error('Inputs has at least one NaN value!')
        end

        if sum(isnan(targets(:))) > 0
            error('Targets has at least one NaN value!')
        end
        
        % Net properties
        x = inputs';
        t = targets';       
        trainFcn = 'trainbr';
        hiddenLayerSize = 1;
        performFcn = 'mse';
        numLayers = numel(hiddenLayerSize);
        
        % Net object
        net = patternnet(hiddenLayerSize, trainFcn, performFcn);
        
        % Layers structure        
        for i = 1:numLayers
            net.layers{i}.transferFcn = 'logsig';
        end
        net.layers{numLayers+1}.transferFcn = 'logsig';
        
        % Setup Division of Data for Training, Validation, Testing
        net.divideFcn = 'dividerand';  % Divide data randomly
        net.divideMode = 'sample';  % Divide up every sample
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;

        % Number of epochs
        net.trainParam.epochs = 1000;
        
        % Training
        [net, tr] = train(net,x,t);
        
        % View the Network
        view(net)
        
        % Test the Network
        y = net(x);
        e = gsubtract(t,y);
        performance = perform(net,t,y);
        tind = vec2ind(t);
        yind = vec2ind(y);
        percentErrors = sum(tind ~= yind)/numel(tind);
        
        % Training Confusion Plot Variables
        yTrn = net(x(:,tr.trainInd));
        tTrn = t(:,tr.trainInd);

        % Validation Confusion Plot Variables
        yVal = net(x(:,tr.valInd));
        tVal = t(:,tr.valInd);

        % Test Confusion Plot Variables
        yTst = net(x(:,tr.testInd));
        tTst = t(:,tr.testInd);

        % Overall Confusion Plot Variables
        yAll = net(x);
        tAll = t;

        % Plot Confusion
        figure, plotconfusion(tTrn, yTrn, 'Training', ...
                              tVal, yVal, 'Validation', ...
                              tTst, yTst, 'Test', ...
                              tAll, yAll, 'Overall')
                          
        wb = formwb(net,net.b,net.iw,net.lw);
        [B,IW,LW] = separatewb(net,wb);
        IW = IW{1,1};
    catch ME
		errorMessage = sprintf('Error in GetWeights().\n The error reported by MATLAB is:\n\n%s', ME.message);
		uiwait(warndlg(errorMessage));
		set(handles.txtInfo, 'String', errorMessage);
    end
	return;
