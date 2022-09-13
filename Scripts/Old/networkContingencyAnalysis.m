%%                           Loading data
%Group 1 data
%load('D:\NetworkConnectivity\RestData\ADHD\ADHDNetwork_NoNone_allSubjStack.mat');
load('D:\NetworkConnectivity\TicData\Tic\TicNetwork_allSubjStack.mat');
allConnectivityStack1   = allConnectivityStack;
finallySelectedEdgeIdx1 = finallySelectedEdgeIdx;

%Group 2 data
%load('D:\NetworkConnectivity\RestData\Control\ControlNetwork_NoNone_allSubjStack.mat');
load('D:\NetworkConnectivity\TicData\Rest\RestNetwork_allSubjStack.mat');
allConnectivityStack2   = allConnectivityStack;
finallySelectedEdgeIdx2 = finallySelectedEdgeIdx;

%Load in network/parcel labels
%load('D:\NetworkConnectivity\RestData\Control\ControlNetwork_NoNone_dipolePairDensity.mat','parcelInfo');
load('D:\NetworkConnectivity\TicData\Tic\TicNetwork_dipolePairDensity.mat','parcelInfo');

%%                  Reorganization of matrix blocks/plotting
%Reorder/group connections based on networks
networkNames = unique(parcelInfo(:,6));

%Exclude any networks (e.g., excludeNetworks = [5])
excludeNetworks = [];
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

%%                         Significance Testing
% Get connections which have sufficient subjects within both groups
overlapEdge = intersect(finallySelectedEdgeIdx1,finallySelectedEdgeIdx2);

%Run significance test on each connection between groups
connGroup1Reshape = reshape(allConnectivityStack1,[size(allConnectivityStack1,1)*size(allConnectivityStack1,2),size(allConnectivityStack1,5)])';
connGroup2Reshape = reshape(allConnectivityStack2,[size(allConnectivityStack2,1)*size(allConnectivityStack2,2),size(allConnectivityStack2,5)])';
connGroup1Reshape(connGroup1Reshape==0)= NaN;
connGroup2Reshape(connGroup2Reshape==0)= NaN;
[~, pValues, ~, stats] = ttest2(connGroup1Reshape, connGroup2Reshape);

pMatrix = reshape(pValues,[size(allConnectivityStack1,1),size(allConnectivityStack1,2)]);
tMatrix = reshape(stats.tstat,[size(allConnectivityStack1,1),size(allConnectivityStack1,2)]);

%Get mask and apply to tstats
pMatrixThresh = pMatrix < 0.05;
tMatrixThresh = tMatrix.*pMatrixThresh;
tMatrixThresh(isnan(tMatrixThresh)) = 0;

%Signed matrix for binary greater than/less than connectivity b/w groups
tMatrixThreshSign = sign(tMatrixThresh); 


%%                      Permutations on each cell
%
bothGroupConnReshape = [connGroup1Reshape;connGroup2Reshape];
numTests = 1000;
nullMatrix = zeros(size(networkNames,1),size(networkNames,1),numTests);

for i = 1:numTests
    permIdx = randperm(size(connGroup1Reshape,1) + size(connGroup2Reshape,1));
    
    permGroup1Conn = bothGroupConnReshape(permIdx(1:size(connGroup1Reshape,1)),:);
    permGroup2Conn = bothGroupConnReshape(permIdx(size(connGroup1Reshape,1)+1:end),:);
    
    [~, permpValues, ~, permstats] = ttest2(permGroup1Conn, permGroup2Conn);
    
    %Use division line idx to calculate # of significant edges
    permpMatrix = reshape(permpValues,[size(allConnectivityStack1,1),size(allConnectivityStack1,2)]);
    permpMatrixThresh = permpMatrix < 0.05;
    
    %Create network # by network # by # iterations matrix for 13*13 null
    %distributions
    for row = 1:length(divisionLineIdx)-1
        parcelBounds1 = [divisionLineIdx(row)+1:divisionLineIdx(row+1)];
        for col = 1:length(divisionLineIdx)-1
            parcelBounds2 = [divisionLineIdx(col)+1:divisionLineIdx(col+1)];
            nullMatrix(row,col,i) = sum(sum(permpMatrixThresh(parcelBounds1,parcelBounds2)));
        end
    end

    if mod(i,5) == 0
        fprintf('Finished permutation iteration %d\n',i)
    end
end
%}


%%                              Plotting

%Plot surviving connections
figure; imagesc(tMatrixThreshSign); set(gcf,'color','w')

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

%Find how many signif connections are present within/between each network
trueMatrix = zeros(size(networkNames,1),size(networkNames,1));
signCount = zeros(size(networkNames,1),size(networkNames,1),2);
for row = 1:length(divisionLineIdx)-1
    parcelBounds1 = [divisionLineIdx(row)+1:divisionLineIdx(row+1)];
    for col = 1:length(divisionLineIdx)-1
        parcelBounds2 = [divisionLineIdx(col)+1:divisionLineIdx(col+1)];
        trueMatrix(row,col) = sum(sum(pMatrixThresh(parcelBounds1,parcelBounds2)));
        
        %How many connections favor each group
        signCount(row,col,1) = length(find(tMatrixThreshSign(parcelBounds1,parcelBounds2)==1));  %Higher conn in group1
        signCount(row,col,2) = length(find(tMatrixThreshSign(parcelBounds1,parcelBounds2)==-1)); %Higher conn in group2
    end
end

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

%Coloring significant (surviving) networks of matrix
for i = 1:size(networkIdx,1)
    hold on
    group1Num = signCount(networkIdx(i,1),networkIdx(i,2),1);
    group2Num = signCount(networkIdx(i,1),networkIdx(i,2),2);
    colorFactor = tanh(2*(group1Num - group2Num)/(group1Num + group2Num));
    boxColor = [0.5 + 0.5*colorFactor , 0, 0.5 - 0.5*colorFactor];
    patch([divisionLineIdx(networkIdx(i,2)), divisionLineIdx(networkIdx(i,2)), divisionLineIdx(networkIdx(i,2)+1), divisionLineIdx(networkIdx(i,2)+1)],... %x coordinates (box)
          [divisionLineIdx(networkIdx(i,1)), divisionLineIdx(networkIdx(i,1)+1), divisionLineIdx(networkIdx(i,1)+1), divisionLineIdx(networkIdx(i,1))],... %y coordinates (box)
          boxColor, 'FaceAlpha', 0.4)
end    