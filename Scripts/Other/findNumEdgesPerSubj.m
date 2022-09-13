function [subjByNetwork1, subjByNetwork2] = findNumEdgesPerSubj(allConnectivityStack1, allConnectivityStack2, divisionLineIdx)
    avgData1 = squeeze(mean(allConnectivityStack1,4));
    subjByNetwork1 = zeros(size(divisionLineIdx,2)-1, size(divisionLineIdx,2)-1, size(avgData1,3));

    for row = 1:length(divisionLineIdx)-1
        parcelBounds1 = [divisionLineIdx(row)+1:divisionLineIdx(row+1)];
        for col = 1:length(divisionLineIdx)-1
            parcelBounds2 = [divisionLineIdx(col)+1:divisionLineIdx(col+1)];
            subjByNetwork1(row,col,:) = sum(sum(avgData1(parcelBounds1,parcelBounds2,:),1),2);
        end
    end
    
    
    avgData2 = squeeze(mean(allConnectivityStack2,4));
    subjByNetwork2 = zeros(size(divisionLineIdx,2)-1, size(divisionLineIdx,2)-1, size(avgData2,3));
    
    for row = 1:length(divisionLineIdx)-1
        parcelBounds1 = [divisionLineIdx(row)+1:divisionLineIdx(row+1)];
        for col = 1:length(divisionLineIdx)-1
            parcelBounds2 = [divisionLineIdx(col)+1:divisionLineIdx(col+1)];
            subjByNetwork2(row,col,:) = sum(sum(avgData2(parcelBounds1,parcelBounds2,:),1),2);
        end
    end
    
    figure; count=1;
    for r = 1:6
        for c = 1:6
            subplot(6,6,count); bar(squeeze(subjByNetwork2(r,c,:))); hold on; bar(squeeze(subjByNetwork1(r,c,:)),'r');
            count = count + 1;
        end
    end
    
    subjByNetwork1 = reshape(subjByNetwork1,[size(subjByNetwork1,1)*size(subjByNetwork1,2),size(subjByNetwork1,3)])';
    subjByNetwork2 = reshape(subjByNetwork2,[size(subjByNetwork2,1)*size(subjByNetwork2,2),size(subjByNetwork2,3)])';
    
    [~, pValues, ~, stats] = ttest2(subjByNetwork1, subjByNetwork2);
    
    pValuesReshape = reshape(pValues,[size(divisionLineIdx,2)-1,size(divisionLineIdx,2)-1]);
    1