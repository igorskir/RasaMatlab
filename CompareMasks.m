function [ acc ] = CompareMasks(mask1,mask2)
%returns Jaccard and Dice assesments of two masks
%mask1 - result of segmentation
%mask2 - original mask

intersectArea=0;
unionArea=0;
for i=1:size(mask1,1)
    for j=1:size(mask1,2)
        if(mask1(i,j)~=0 || mask2(i,j)~=0)
            unionArea=unionArea+1;
        end
        if(mask1(i,j)~=0 && mask2(i,j)~=0)
            intersectArea=intersectArea+1;
        end
    end
end
acc.jac=intersectArea/unionArea;
acc.dice=2*acc.jac/(1+acc.jac);
end

