function res = DensitiesDifference(cathDataArr, tissueDataArr, nBins)
%HOW TO USE
%����� �������� workspace � data.mat, ��� ���� ������ �� ����, ����������
%��� �������, ��� ������ �� ������ � ��������
%� ������� ����� �������:
%r = DensitiesDifference(dataCatheter(:,2), dataTissue(:,2),300);
%������� ���� �� ���� ������ Nx1, � �� �������. 

% [cathBinElements,edges] = histcounts(cathMatrix(:,2), 'BinMethod', 'fd');
% [tissBinElements,edges] = histcounts(tissueMatrix(:,2), 'BinMethod', 'fd');
    
    if size(cathDataArr,2) > 1
        cathDataArr = cathDataArr';
    end
    if size(tissueDataArr,2) > 1
        tissueDataArr = tissueDataArr'; 
    end
    
    %[cathBinElements,edges] = histcounts(cathDataArr,nBins);
    %[tissBinElements,edges] = histcounts(tissueDataArr, nBins);
    
    
    
    
    cathBinElementsNorm = cathBinElements/size(cathDataArr,1);
    tissBinElementsNorm = tissBinElements/size(tissueDataArr,1);
    
    res = 0;
    for i = 1:nBins
        currentDelta = cathBinElementsNorm(i) - tissBinElementsNorm(i);
        if cathBinElements(i) ~= 0 && currentDelta > 0
            res = res + currentDelta; 
        end    
    end
    
    %����������� ������������    
    plot(cathBinElementsNorm);
    hold on;
    plot(tissBinElementsNorm);
    hold off;
    legend('Catheter PDF','Tissue PDF')

end