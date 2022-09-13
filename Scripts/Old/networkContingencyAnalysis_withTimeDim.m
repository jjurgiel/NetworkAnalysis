%
%%                           Loading data
%Group 1 data
%load('D:\NetworkConnectivity\RestData\ADHD\ADHDNetwork_NoNone_allSubjStack.mat');
%load('D:\NetworkConnectivity\TicData\Rest\Rest30p_allSubjStack.mat');
load('D:\NetworkConnectivity\FlankerData\Control\Incongruent\IncongruentNetwork30p_allSubjStack.mat');
allConnectivityStack1   = allConnectivityStack;
finallySelectedEdgeIdx1 = finallySelectedEdgeIdx;

%Group 2 data
%load('D:\NetworkConnectivity\RestData\Control\ControlNetwork_NoNone_allSubjStack.mat');
%load('D:\NetworkConnectivity\TicData\Tic\Tic30p_allSubjStack.mat');
load('D:\NetworkConnectivity\FlankerData\Tic\Incongruent\IncongruentNetwork30p_allSubjStack.mat');
allConnectivityStack2   = allConnectivityStack;
finallySelectedEdgeIdx2 = finallySelectedEdgeIdx;

%Load in network/parcel labels
%load('D:\NetworkConnectivity\RestData\Control\ControlNetwork_NoNone_dipolePairDensity.mat','parcelInfo');
%load('D:\NetworkConnectivity\TicData\Tic\Tic30p_dipolePairDensity.mat','parcelInfo');
load('D:\NetworkConnectivity\FlankerData\Tic\Incongruent\IncongruentNetwork30p_dipolePairDensity.mat','parcelInfo');


%Other useful info

%baseline = [-3.5 -2.5];
baseline = [-0.55 -0.05];
numTests = 500;

%%                  Reorganization of matrix blocks/plotting
%Reorder/group connections based on networks
networkNames = unique(parcelInfo(:,6));

%Exclude any networks (e.g., excludeNetworks = [5])
excludeNetworks               = [];
networkNames(excludeNetworks) = [];

newParcelOrder = []; divisionLineIdx = 0; 
for i = 1:size(networkNames,1)
    parcelsIdxNetwork = find(strcmp(parcelInfo(:,6),networkNames(i))==1);
    newParcelOrder = [newParcelOrder ; parcelsIdxNetwork];
    
    %Division lines for plotting later
    divisionLineIdx(i+1) = divisionLineIdx(i) + length(parcelsIdxNetwork);
end


allConnectivityStack1 = allConnectivityStack1(newParcelOrder',newParcelOrder',:,:,:);
allConnectivityStack2 = allConnectivityStack2(newParcelOrder',newParcelOrder',:,:,:);

%Other things
%[subjByNetwork1, subjByNetwork2] = findNumEdgesPerSubj(logical(allConnectivityStack1), logical(allConnectivityStack2), divisionLineIdx);

%# subjects in each network
%Adjust finallySelectedEdge (because we re-ordered matrix)
edgeMatrix1 = zeros(size(allConnectivityStack1,1),size(allConnectivityStack1,2));
edgeMatrix1(finallySelectedEdgeIdx1) = 1;
finallySelectedEdgeIdx1New = find(edgeMatrix1(newParcelOrder',newParcelOrder')>0);

edgeMatrix2 = zeros(size(allConnectivityStack2,1),size(allConnectivityStack2,2));
edgeMatrix2(finallySelectedEdgeIdx2) = 1;
finallySelectedEdgeIdx2New = find(edgeMatrix2(newParcelOrder',newParcelOrder')>0);

%%                         Baseline data (subtraction)
baselineIdx = find(latencies > baseline(1) & latencies < baseline(2));
baseline1   = mean(allConnectivityStack1(:,:,1,baselineIdx,:),4);
allConnectivityStack1 = allConnectivityStack1 - baseline1;

baseline2   = mean(allConnectivityStack2(:,:,1,baselineIdx,:),4);
allConnectivityStack2 = allConnectivityStack2 - baseline2;
disp('Connectivity data baselined')

%%                  Array setup for significance Testing
% Get connections which have sufficient subjects within both groups
overlapEdge              = intersect(finallySelectedEdgeIdx1New,finallySelectedEdgeIdx2New);

%Set up true and null arrays in memory
tMatrixThreshSign        = int8(zeros(size(allConnectivityStack1,1),size(allConnectivityStack1,2),size(allConnectivityStack1,4)));
networkCount             = single(zeros(size(networkNames,1),size(networkNames,1),size(allConnectivityStack1,4)));
networkPercent           = single(networkCount);

nullNetworkCount         = single(zeros(size(networkNames,1),size(networkNames,1),size(allConnectivityStack1,4),numTests));
nullNetworkPercent       = single(nullNetworkCount);


%%                         Significance Testing
connGroup1Reshape = squeeze(reshape(allConnectivityStack1,[size(allConnectivityStack1,1)*size(allConnectivityStack1,2),size(allConnectivityStack1,3),size(allConnectivityStack1,4),size(allConnectivityStack1,5)]));
connGroup2Reshape = squeeze(reshape(allConnectivityStack2,[size(allConnectivityStack2,1)*size(allConnectivityStack2,2),size(allConnectivityStack2,3),size(allConnectivityStack2,4),size(allConnectivityStack2,5)]));
connGroup1Reshape = connGroup1Reshape(overlapEdge,:,:);
connGroup2Reshape = connGroup2Reshape(overlapEdge,:,:);
connGroup1Reshape(connGroup1Reshape==0) = NaN;
connGroup2Reshape(connGroup2Reshape==0) = NaN;
connGroup1Reshape = reshape(connGroup1Reshape,[size(connGroup1Reshape,1)*size(connGroup1Reshape,2),size(connGroup1Reshape,3)])';
connGroup2Reshape = reshape(connGroup2Reshape,[size(connGroup2Reshape,1)*size(connGroup2Reshape,2),size(connGroup2Reshape,3)])';

pValuesAll = single(ones(size(tMatrixThreshSign)));
pValuesAll = reshape(pValuesAll,[size(pValuesAll,1)*size(pValuesAll,2),size(pValuesAll,3)]);
tMatrixAll = pValuesAll;

%Run t-test
[~, pValues, ~, stats] = ttest2(connGroup1Reshape, connGroup2Reshape);

pValuesReshape = reshape(pValues,[size(overlapEdge,1),size(allConnectivityStack2,4)]);
pValuesAll(overlapEdge,:) = pValuesReshape;
pValuesAll = reshape(pValuesAll,size(tMatrixThreshSign)); 

tMatrixReshape = reshape(stats.tstat,[size(overlapEdge,1),size(allConnectivityStack2,4)]);
tMatrixAll(overlapEdge,:) = tMatrixReshape;
tMatrixAll = reshape(tMatrixAll,size(tMatrixThreshSign)); 

pValuesAll = pValuesAll < 0.05;
tMatrixAllThresh = tMatrixAll .* pValuesAll;
tMatrixAllThresh(isnan(tMatrixAllThresh)) = 0;
tMatrixThreshSign = sign(tMatrixAllThresh); 


%% Calculate details about each network pair (available conn, counts, etc)
idxMatrix     = zeros(size(parcelInfo,1),size(parcelInfo,1));
edgeAvailable = zeros(size(networkNames,1),size(networkNames,1));
val = 1;
for row = 1:length(divisionLineIdx)-1
    parcelBounds1 = [divisionLineIdx(row)+1:divisionLineIdx(row+1)];
    for col = 1:length(divisionLineIdx)-1
        parcelBounds2 = [divisionLineIdx(col)+1:divisionLineIdx(col+1)];
        
        %Set up matrix of network idx
        idxMatrix(parcelBounds1, parcelBounds2) = val;
        
        %Find # available edges in each network
        networkIdx = find(idxMatrix == val);
        edgeAvailable(val) = length(intersect(overlapEdge,networkIdx));
        
        %Calculate network counts and percents
        networkMatrix = abs(tMatrixThreshSign(parcelBounds1,parcelBounds2,:));
        networkCount(row,col,:) = sum(sum(networkMatrix,1),2);
        networkPercent(row,col,:) = networkCount(row,col,:)/edgeAvailable(val);
        
        val = val + 1;
    end
end

fprintf('\nFinished True Statistics\n')

%%                      Permutations on each cell
%Group data
bothGroupConnReshape = [connGroup1Reshape;connGroup2Reshape];

for t = 1:numTests
    %Permute data
    permIdx        = randperm(size(connGroup1Reshape,1) + size(connGroup2Reshape,1));
    permGroup1Conn = bothGroupConnReshape(permIdx(1:size(connGroup1Reshape,1)),:);
    permGroup2Conn = bothGroupConnReshape(permIdx(size(connGroup1Reshape,1)+1:end),:);

    [~, permpValues, ~, permstats] = ttest2(permGroup1Conn, permGroup2Conn);

    permpValuesReshape = reshape(permpValues,[size(overlapEdge,1),size(allConnectivityStack2,4)]);
    %Initiate
    permpValuesAll = single(ones(size(tMatrixThreshSign,1)*size(tMatrixThreshSign,2),size(tMatrixThreshSign,3)));
    permpValuesAll(overlapEdge,:) = permpValuesReshape;
    permpValuesAll = reshape(permpValuesAll,size(tMatrixThreshSign)); 
    
    %Threshold
    permpValuesAll = permpValuesAll < 0.05;
    
    %Create network x network x time x iteration perm distribution
    for row = 1:length(divisionLineIdx)-1
        parcelBounds1 = [divisionLineIdx(row)+1:divisionLineIdx(row+1)];
        for col = 1:length(divisionLineIdx)-1
            parcelBounds2 = [divisionLineIdx(col)+1:divisionLineIdx(col+1)];
            networkMatrix = abs(permpValuesAll(parcelBounds1,parcelBounds2,:));
            nullNetworkCount(row,col,:,t)   = sum(sum(networkMatrix,1),2);
            nullNetworkPercent(row,col,:,t) = nullNetworkCount(row,col,:,t)/edgeAvailable(row, col);
        end
    end

    if mod(t,5) == 0
        fprintf('Finished permutation iteration %d of %d \n', t, numTests)
    end

    fprintf('.')
end
fprintf('\nFinished perm statistics\n')
%}

%Grab 95th percentile for each network (for both count and percent)
nullNetworkCountReshape   = reshape(nullNetworkCount,size(nullNetworkCount,1),size(nullNetworkCount,2),size(nullNetworkCount,3)*size(nullNetworkCount,4));
nullNetworkPercentReshape = reshape(nullNetworkPercent,size(nullNetworkCount,1),size(nullNetworkCount,2),size(nullNetworkCount,3)*size(nullNetworkCount,4));
countNetworkThresholds = prctile(nullNetworkCountReshape,95,3); %95% threshold over 3rd dim
percentNetworkThresholds = prctile(nullNetworkPercentReshape,95,3); %95% threshold over 3rd dim
percentNetworkThresholdGlobal = prctile(nullNetworkPercentReshape(:),95); %95% threshold over 3rd dim

%For JD: Plotting network % threshold as a function of edges available
figure; scatter(percentNetworkThresholds(:),edgeAvailable(:))
hold on; xline(percentNetworkThresholdGlobal,'-','Global % Thresh')
ylabel('# Edges Available in given network')
xlabel('95th percentile of % availalbe edges significant null distribution')

%Apply masks to true data 
networkCountMasked = networkCount > countNetworkThresholds;
networkPercentMasked = networkPercent > percentNetworkThresholds;
networkPercentGlobalMasked = networkPercent > percentNetworkThresholdGlobal;

%Apply masks to permutation data
nullNetworkCountMasked = nullNetworkCount > countNetworkThresholds;
nullNetworkPercentMasked = nullNetworkPercent > percentNetworkThresholds;
nullNetworkPercentGlobalMasked = nullNetworkPercent > percentNetworkThresholdGlobal;

%Now create null of "count/percent cluster masses"
%This is for figuring out how many sequentual time bines is needed for
%significance
approach = 'percentLocal'; %methods: percentGlobalTime, percentGlobal, percentLocal, countLocal
% percentGlobalTime: sums up how many adjacent time points
% percentGlobal: finds max sum percent across all networks
% percentLocal: finds max sum percent for each network (needs FDR)
% countLocal: finds max sum count for each network (needs FDR)
permArray = single(zeros(size(nullNetworkCount,1),size(nullNetworkCount,2),size(nullNetworkCount,4)));
for i = 1:size(nullNetworkCount,4)
    for row = 1:length(divisionLineIdx)-1
        for col = 1:length(divisionLineIdx)-1
            switch approach
                case 'percentGlobalTime'
                    iterData =  squeeze(nullNetworkPercentGlobalMasked(row,col,:,i));
                    nullData = squeeze(nullNetworkPercent(row,col,:,i));
                case 'percentGlobal'
                    iterData =  squeeze(nullNetworkPercentGlobalMasked(row,col,:,i));
                    nullData = squeeze(nullNetworkPercent(row,col,:,i));
                case 'percentLocal'
                    iterData = squeeze(nullNetworkPercentMasked(row,col,:,i));
                    nullData = squeeze(nullNetworkPercent(row,col,:,i));
                case 'countLocal'
                    iterData = squeeze(nullNetworkCountMasked(row,col,:,i));
                    nullData = squeeze(nullNetworkCount(row,col,:,i));                    
            end
                [labeledData, n] = bwlabel(iterData);
                if n > 0
                    sumOfBlobs = [];
                    for label = 1:n
                        switch approach
                            case 'percentGlobalTime'
                                sumOfBlobs(label) = length(find(labeledData==label));
                            case {'percentGlobal', 'percentLocal', 'countLocal'}
                                sumOfBlobs(label) = sum(nullData(labeledData==label));
                        end                   
                    end
                    permArray(row,col,i) = max(sumOfBlobs);
                    1;
                end
        end
    end
end

%Get permutation-estimated threshold
switch approach
    case {'percentGlobalTime', 'percentGlobal'}
        permArrayReshape = reshape(permArray,[size(permArray,1)*size(permArray,2),size(permArray,3)]);
        permArrayReshapeSorted = sort(permArrayReshape,2,'descend');
        permThresh = prctile(permArrayReshapeSorted(:,1),95); %Time threshold of null
        permThresh = repmat(permThresh,size(permArray,1),size(permArray,2)); 
    case {'percentLocal', 'countLocal'}
        permThresh = prctile(permArray,95,3);
end

maxBlob = 0;
%Now apply time threshold to true data
for row = 1:length(divisionLineIdx)-1
    for col = 1:length(divisionLineIdx)-1
            switch approach
                case 'percentGlobalTime'
                    iterData =  squeeze(networkPercentGlobalMasked(row,col,:));
                    trueData = squeeze(networkPercent(row,col,:));
                case 'percentGlobal'
                    iterData =  squeeze(networkPercentGlobalMasked(row,col,:));
                    trueData = squeeze(networkPercent(row,col,:));
                case 'percentLocal'
                    iterData = squeeze(networkPercentMasked(row,col,:));
                    trueData = squeeze(networkPercent(row,col,:));
                case 'countLocal'
                    iterData = squeeze(networkCountMasked(row,col,:));
                    trueData = squeeze(networkCount(row,col,:));                    
            end
            [labeledData, n] = bwlabel(iterData);
            if n > 0
                for label = 1:n
                    switch approach
                        case 'percentGlobalTime'
                            sumOfBlobs = length(find(labeledData==label));
                        case {'percentGlobal', 'percentLocal', 'countLocal'}
                            sumOfBlobs = sum(trueData(labeledData==label));
                    end    
                    maxBlob = max(maxBlob, sumOfBlobs);
                    if sumOfBlobs < permThresh(row,col)
                        iterData((labeledData==label)) = 0;
                    else
                        1;
                    end
                end
                networkMaskedCorrected(row,col,:) = iterData;
            end
    end
end
%}
if ~any(networkMaskedCorrected(:))
    error('No surviving Results')
else
%%              Plotting network x time logical array 
    networkMaskedCorrectedReshape = reshape(networkMaskedCorrected,size(networkMaskedCorrected,1)*size(networkMaskedCorrected,2),size(networkMaskedCorrected,3));
    [rowIdx, colIdx] = find(any(networkMaskedCorrected,3));
    %[rowIdx, colIdx] = find(sum(networkMaskedCorrected,3)>=0);
    connLabels = {}; plotArray = [];
    for i = 1:length(rowIdx)
        connLabels{i} = horzcat(networkNames{colIdx(i)},' to ',networkNames{rowIdx(i)});
        plotArray = [plotArray;networkMaskedCorrected(rowIdx(i),colIdx(i),:)];
    end
    plotArray = squeeze(plotArray);
    %Plotting/axes fixes
    figure; imagesc(squeeze(plotArray))
    yTicks = [1:size(connLabels,2)];
    set(gca, 'YTick', yTicks, 'YTickLabel', connLabels)

    xTickLabels = latencies(1):.25:latencies(end)+0.25;
    xTicks = linspace(1, size(plotArray,2), numel(xTickLabels));
    set(gca, 'XTick', xTicks, 'XTickLabel', xTickLabels)
end
%%                              Plotting

%FOR JD 5/17: 
%Permutation distribution of local percents
figure
for i = 1:13
    subplot(4,4,i)
    plotData = squeeze(nullNetworkPercent(i,i,:,:));
    hist(plotData(:),100);
    title(sprintf('N: %s, Avail: %d, 95p: %s',networkNames{i},edgeAvailable(i,i), num2str(prctile(plotData(:),95))))
end
%
plotArray = squeeze(plotArray); midValue = [];
signifIdx = {};
for i = 1:length(rowIdx)
    [labeledIdx, n] = bwlabel(plotArray(i,:));
    for j = 1:n
        signifIdx{i} = find(labeledIdx == j);
        midValue = [midValue; [i round(median(signifIdx{i}))]];
    end
     
end
   

%Plotting significant results (before FDR correction)
trueMatrix = zeros(size(networkNames,1),size(networkNames,1));
edgeFull = zeros(size(allConnectivityStack1,1),size(allConnectivityStack1,2),size(midValue,1));
signCount = zeros(size(midValue));
for i = 1:size(midValue,1)
    parcelBounds1 = [divisionLineIdx(rowIdx(i))+1:divisionLineIdx(rowIdx(i)+1)];
    parcelBounds2 = [divisionLineIdx(colIdx(i))+1:divisionLineIdx(colIdx(i)+1)];
    
    signCount(i,1) = length(find(tMatrixThreshSign(parcelBounds1,parcelBounds2,midValue(i,2))==1));  %Higher conn in group1
    signCount(i,2) = length(find(tMatrixThreshSign(parcelBounds1,parcelBounds2,midValue(i,2))==-1)); %Higher conn in group2    
    
    
    %Plot surviving connections
    figure; imagesc(tMatrixThreshSign(:,:,midValue(i,2))); set(gcf,'color','w')

    %Draw network division lines
    arrayfun(@(a)xline(a),divisionLineIdx(2:end-1))    
    arrayfun(@(a)yline(a),divisionLineIdx(2:end-1))

    %Plot diagonal
    x = [1, size(allConnectivityStack1,2)]; y = [1, size(allConnectivityStack1,1)];
    line(x,y,'color',[0 0 0])

    %Get indices for labeling plot (middle of network coordinates)
    labelCoords = round(mean([divisionLineIdx(1:end-1);divisionLineIdx(2:end)]));
    xticks(labelCoords); xticklabels(networkNames); xtickangle(45);
    yticks(labelCoords); yticklabels(networkNames); 
    set(gca,'TickLength',[0 0]);

    %Color shade connections based on direction of group difference
    colormap([0 0 1 ; 1 1 1; 1 0 0])
    set(gca,'XAxisLocation','top','YAxisLocation','left');

    %Make matrix square (visualization purpose)
    x0=100; y0=100; width=1000; height=800; set(gcf,'position',[x0,y0,width,height]); axis square
    xlabel('From Region'); ylabel('To Region')
    
    %Shading square
    group1Num = signCount(i,1);
    group2Num = signCount(i,2);
    colorFactor = tanh(2*(group1Num - group2Num)/(group1Num + group2Num));
    
    boxColor = [0.5 + 0.5*colorFactor , 0, 0.5 - 0.5*colorFactor];
    patch([divisionLineIdx(colIdx(i)), divisionLineIdx(colIdx(i)), divisionLineIdx(colIdx(i)+1), divisionLineIdx(colIdx(i)+1)],... %x coordinates (box)
          [divisionLineIdx(rowIdx(i)), divisionLineIdx(rowIdx(i)+1), divisionLineIdx(rowIdx(i)+1), divisionLineIdx(rowIdx(i))],... %y coordinates (box)
          boxColor, 'FaceAlpha', 0.4)    
      
    annotation('textbox', [0.45, 0.01, 0.1, 0.1], 'string', horzcat('Extracted from time: ',num2str(latencies(min(signifIdx{i}))),' to ',num2str(latencies(max(signifIdx{i}))),'s'))
    
    %Create Edge matrix for connectivity plot
    edgeFull(divisionLineIdx(rowIdx(i))+1:divisionLineIdx(rowIdx(i)+1),divisionLineIdx(colIdx(i))+1:divisionLineIdx(colIdx(i)+1),i) = tMatrixThreshSign(divisionLineIdx(rowIdx(i))+1:divisionLineIdx(rowIdx(i)+1),divisionLineIdx(colIdx(i))+1:divisionLineIdx(colIdx(i)+1),midValue(i,2));
    edgeFull(edgeFull==-1) = 0.5;
    
    

end

1
% signCount = zeros(size(networkNames,1),size(networkNames,1),2);
% for row = 1:length(divisionLineIdx)-1
%     parcelBounds1 = 1;
%     for col = 1:length(divisionLineIdx)-1
%         parcelBounds2 = [divisionLineIdx(col)+1:divisionLineIdx(col+1)];
%         %trueMatrix(row,col) = sum(sum(pMatrixThresh(parcelBounds1,parcelBounds2)));
%         
%         %How many connections favor each group
%         signCount(row,col,1) = length(find(tMatrixThreshSign(parcelBounds1,parcelBounds2,midValue(i))==1));  %Higher conn in group1
%         signCount(row,col,2) = length(find(tMatrixThreshSign(parcelBounds1,parcelBounds2,midValue(i))==-1)); %Higher conn in group2
%     end
% end

1;

%FDR Correction on survivingNetworks
networkThresholds = prctile(nullMatrix,95,3); %95% threshold over 3rd dim
survivingNetworks = trueMatrix > networkThresholds;
networkIdx = [];
[networkIdx(:,1),networkIdx(:,2)] = find(survivingNetworks==1);
%Correction done below
tableForFDR = []; count = 1;
for i = 1:size(networkNames,1)
    for j = 1:size(networkNames,1)
        nullForNetwork = squeeze(nullMatrix(i,j,:));
        totalSignifConn = sum(signCount(i,j,:));
        pValueForNetwork = sum(nullForNetwork>totalSignifConn)/numTests;
        
        tableForFDR(count,1) = i;
        tableForFDR(count,2) = j;
        tableForFDR(count,3) = pValueForNetwork;
        
        count = count + 1;
    end
end
%Get adjusted p values
tableForFDRSorted = sortrows(tableForFDR,3);
for i = 1:length(tableForFDRSorted)
    tableForFDRSorted(i,4) = tableForFDRSorted(i,3)*length(tableForFDRSorted)/(i);
end




