%
%%                      Selecting files to Load
%{
%Select group 1 data
[group1Files, group1Folder] = uigetfile('*.conn', 'MultiSelect', 'on');
if ~any(group1Folder); disp('Cancelled.'); return; end

%Select group 2 data
[group2Files, group2Folder] = uigetfile('*.conn', 'MultiSelect', 'on');
if ~any(group1Folder); disp('Cancelled.'); return; end

%
%%                          Network Analysis

%Will run analysis network-by-network, in order to allow for large
%datasets that cant store entire matrix for all subj in memory

%First get network indices
load(horzcat(group1Folder, group1Files{1}), '-mat', 'parcelInfo', 'frequencies', 'latencies');

networkNames = unique(parcelInfo(:, 6));

%Exclude any networks (e.g., excludeNetworks = [5])
excludeNetworks               = [];
networkNames(excludeNetworks) = [];

parcelsIdxNetwork = {}; divisionLineIdx = 0; 
totalParcels = 0;
for i = 1:size(networkNames,1)
    parcelsIdxNetwork{i} = find(strcmp(parcelInfo(:,6), networkNames(i))==1);
    totalParcels = totalParcels + length(parcelsIdxNetwork{i});
    %Division lines for plotting later
    divisionLineIdx(i+1) = divisionLineIdx(i) + length(parcelsIdxNetwork{i});
end
%

% Now go through all network pairs and perform analysis on each
numTests = 200;
countByNetwork = int8(zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), length(frequencies), length(latencies))); %This stores counts of each network for each time/freq bin
countByNetworkPerm = int8(zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), length(frequencies), length(latencies), numTests)); %Store of counts for perm tests

pValuesAll = uint8(ones(totalParcels, totalParcels, length(frequencies), length(latencies))); %This stores all p values for entire matrix (for plotting mainly)
%tStatsSignAll = int8(zeros(totalParcels,totalParcels,length(frequencies), length(latencies))); %This stores signs for all matrix

countByNetworkMaxStatistic = zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), numTests);

for network1 = 7:size(parcelsIdxNetwork, 2)
    for network2 = 7:size(parcelsIdxNetwork, 2)
        %Coordinates of the pair of networks
        toIdx   = parcelsIdxNetwork{network1};
        fromIdx = parcelsIdxNetwork{network2};
        
        %
        %Load in group 1 data and pull out toIdx x fromIdx x t x f matrix
        connMatSegGroup1 = single(zeros(length(toIdx), length(fromIdx), length(frequencies), length(latencies), size(group1Files,2)));
        for n = 1:size(group1Files, 2)
            load(horzcat(group1Folder, group1Files{n}), '-mat');
            
            %This fcn will get the full list of idx of the stored data
            nonzeroIdxFull = reshape(bsxfun(@plus, double(subjectEdgeIdx), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]),[],1);
            %Then will re-create the original large matrix (n x n x freq x time)
            connMatFull = zeros(dims);
            connMatFull(nonzeroIdxFull) = nonzeroConn;
            connMatSegGroup1(:,:,:,:,n) = single(connMatFull(toIdx, fromIdx, :, :));

            disp(horzcat('[Group 1] Loaded file ', num2str(n), ' for network combo: ', num2str(network1), ' ', num2str(network2)))
        end
        finallySelectedEdgeIdx1 = finallySelectedEdgeIdx;
        
        
        %Load in group 1 data and pull out toIdx x fromIdx x t x f matrix
        connMatSegGroup2 = single(zeros(length(toIdx), length(fromIdx), length(frequencies), length(latencies), size(group1Files, 2)));
        for n = 1:size(group2Files, 2)
            load(horzcat(group2Folder,group2Files{n}),'-mat');
            
            %This fcn will get the full list of idx of the stored data
            nonzeroIdxFull = reshape(bsxfun(@plus, double(subjectEdgeIdx), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]), [], 1);
            %Then will re-create the original large matrix (n x n x freq x time)
            connMatFull = zeros(dims);
            connMatFull(nonzeroIdxFull) = nonzeroConn;
            connMatSegGroup2(:,:,:,:,n) = single(connMatFull(toIdx, fromIdx, :, :));
            
            disp(horzcat('[Group 2] Loaded file ', num2str(n), ' for network combo: ', num2str(network1), ' ', num2str(network2)))
        end
        finallySelectedEdgeIdx2 = finallySelectedEdgeIdx;
        %
        
        % Baseline data
        %baseline = [-3.5 -2.5];
        baseline = [-0.55 -0.05];
        baselineIdx = find(latencies > baseline(1) & latencies < baseline(2));
        
        %Group 1 baseline
        baseline1   = mean(connMatSegGroup1(:, :, :, baselineIdx,:),4);
        connMatSegGroup1 = connMatSegGroup1 - baseline1;
        
        %Group 2 baseline
        baseline2   = mean(connMatSegGroup2(:, :, :, baselineIdx,:),4);
        connMatSegGroup2 = connMatSegGroup2 - baseline2;
        
        disp('Connectivity data baselined')

        %%                  Significance Testing
        % Get connections which have sufficient subjects within both groups
        overlapEdge = intersect(finallySelectedEdgeIdx1, finallySelectedEdgeIdx2);  
        overlapEdgeMat = zeros(dims(1), dims(2));
        overlapEdgeMat(overlapEdge) = 1;
        connMatSegGroup1Pair = connMatSegGroup1 .* overlapEdgeMat(toIdx, fromIdx);
        connMatSegGroup2Pair = connMatSegGroup2 .* overlapEdgeMat(toIdx, fromIdx);
        
        %Find non-zero edges to test
        selectedEdges = find(overlapEdgeMat(toIdx, fromIdx) == 1);

        %Select Edges meeting subj criteria for this network
        connMatSegGroup1Reshape = reshape(connMatSegGroup1Pair, [size(connMatSegGroup1Pair,1)*size(connMatSegGroup1Pair,2), size(connMatSegGroup1Pair,3), size(connMatSegGroup1Pair,4), size(connMatSegGroup1Pair,5)]);
        connMatSegGroup2Reshape = reshape(connMatSegGroup2Pair, [size(connMatSegGroup2Pair,1)*size(connMatSegGroup2Pair,2), size(connMatSegGroup2Pair,3), size(connMatSegGroup2Pair,4), size(connMatSegGroup2Pair,5)]);
        
        connMatSegGroup1Selected = connMatSegGroup1Reshape(selectedEdges,:,:,:);
        connMatSegGroup2Selected = connMatSegGroup2Reshape(selectedEdges,:,:,:);
        
        %Statistics on each edge meeting subj Criteria
        [pValuesNetwork, tStatsNetwork] = runStatistics(connMatSegGroup1Selected, connMatSegGroup2Selected);
        
        pValuesNetworkThresh = pValuesNetwork < 0.05;
        pValuesCount = squeeze(sum(pValuesNetworkThresh,1));
        
        %Store pValuesCount in network x network array
        countByNetwork(network1, network2, :, :) = pValuesCount; %**
        
        %Grab direction of significant edges for later use/plotting
       % pValuesTemp = ones(size(connMatSegGroup1Reshape,1), size(connMatSegGroup1Reshape,2), size(connMatSegGroup1Reshape,3));
        tStatsTemp  = zeros(size(connMatSegGroup1Reshape,1), size(connMatSegGroup1Reshape,2), size(connMatSegGroup1Reshape,3));
        
       % pValuesTemp(selectedEdges, :, :) = pValuesNetworkThresh;
        tStatsTemp(selectedEdges, :, :) = tStatsNetwork.*pValuesNetworkThresh;
        
        %pValuesAll(toIdx, fromIdx, :, :) = uint8(reshape(pValuesTemp, [size(connMatSegGroup1Pair,1), size(connMatSegGroup1Pair,2), size(pValuesTemp,2), size(pValuesTemp,3)]));
        tStatsSignAll(toIdx, fromIdx, :, :)  = sign(reshape(tStatsTemp, [size(connMatSegGroup1Pair,1), size(connMatSegGroup1Pair,2), size(tStatsTemp,2), size(tStatsTemp,3)])); %**
        

        %%                      Permutation Testing
        %
        %Put data together to pull out permutation idx
        bothGroupConnReshape = cat(4, connMatSegGroup1Selected, connMatSegGroup2Selected);
        
        pValuesPermCount = zeros(size(bothGroupConnReshape,2), size(bothGroupConnReshape, 3), numTests);
        for t = 1:numTests
            %Permute data
            permIdx        = randperm(size(bothGroupConnReshape,4));
            permGroup1Conn = bothGroupConnReshape(:, :, :, permIdx(1:size(connMatSegGroup1Selected, 4)));
            permGroup2Conn = bothGroupConnReshape(:, :, :, permIdx(size(connMatSegGroup1Selected, 4)+1:end));
            
            %t-tests
            [pValuesNetworkPerm, tStatsNetworkPerm] = runStatistics(permGroup1Conn, permGroup2Conn);
            
            %Mask and get network count 
            pValuesNetworkPermThresh = pValuesNetworkPerm < 0.05;
            pValuesPermCount(:,:,t) = squeeze(sum(pValuesNetworkPermThresh,1));
            
            disp(horzcat('Network combo: ', num2str(network1), ' ', num2str(network2), ' - Finished permutation iteration ', num2str(t), '/', num2str(numTests)))
        end
        
        %Mask perm data using this count (network x network x t x f x perm)
        countByNetworkPerm(network1, network2, :, :, :) = pValuesPermCount; %**
        %
        disp('Finished Network Pair')
    end
end
%
save(['D:\NetworkConnectivity\FlankerData\' 'flanker_tStatistics'],...
    'countByNetwork', 'tStatsSignAll', 'countByNetworkPerm', ...
    'numTests', 'networkNames', 'frequencies', 'latencies', ...
    'parcelsIdxNetwork', 'divisionLineIdx', 'baselineIdx', ...
    'group1Folder', 'group1Files', 'group2Folder', 'group2Files', '-v7.3');
%}

%Save the above to a file for usage on diff fwer control p values...

%{
%%          Now do FWER control and check out significant results
controlledNetworkResults = zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), length(frequencies), length(latencies));
countByNetworkMaxStatistic = zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), size(countByNetworkPerm,5));

%FWER control in each network
for network1 = 1:size(parcelsIdxNetwork, 2)
    for network2 = 1:size(parcelsIdxNetwork, 2)
        
        countMatrix     = single(squeeze(countByNetwork(network1, network2, :, :)));
        countMatrixPerm = single(squeeze(countByNetworkPerm(network1, network2, :, :, :)));
        
        %95th percentile threshold of counts
        countThreshold = prctile(countMatrixPerm(:), 95);
        
        %Mask matrices
        countMatrixMask = countMatrix .* (countMatrix > countThreshold);
        countMatrixPermMask = countMatrixPerm .* (countMatrixPerm > countThreshold);
        
        if any(any(countMatrixMask))
            
            %First get perm distribution of blobs
            for i = 1:size(countMatrixPermMask, 3)
                    countData = countMatrixPermMask(:, :, i);
                    [labeledData, n] = bwlabel(countData);
                    %Find cluster mass of each 'island' of bins
                    if n > 0
                        sumOfBlobs = zeros(1,n);
                        for label = 1:n
                            sumOfBlobs(label) = sum(countData(labeledData == label));                     
                        end
                        %Max statistic
                        countByNetworkMaxStatistic(network1, network2, i) = max(sumOfBlobs);
                    else
                        countByNetworkMaxStatistic(network1, network2, i) = 0;
                    end
            end
            
            %Nth percentile threshold for cluster-level correction
            blobThresh = prctile(countByNetworkMaxStatistic(network1, network2, :), 95);
            
            %Apply blob threshold to True data
            [labeledData, n] = bwlabel(countMatrixMask);
            for label = 1:n
                countSum = sum(countMatrixMask(labeledData == label));    

                if countSum < blobThresh
                    countMatrixMask(labeledData == label) = 0;
                end
            end

            %Store FWER controlled results
            controlledNetworkResults(network1, network2, :, :) = countMatrixMask;
            if any(any(countMatrixMask))
                disp(horzcat('Significant result found in network ', num2str(network1), ' ', num2str(network2)))
            end
        else
            disp('No significant results')
        end
    end
end


%Find regions with significant results
[row, col] = find(any(any(controlledNetworkResults, 3), 4) == 1);

%}

if isempty(row)
    disp('No significant results')
else
    for i = 2:length(row)
        %%          Plotting significant networks & blobs
        resultMatrix = squeeze(controlledNetworkResults(row(i), col(i), :, :));
        
        %Create plot
        edgeFull = plotNCA_Matrix(resultMatrix, row(i), col(i), divisionLineIdx, tStatsSignAll, networkNames, frequencies, latencies);
        
        %%           Extracting single-subject measures
        %Find connections in significant blob
        networkPair = tStatsSignAll(parcelsIdxNetwork{row(i)}, parcelsIdxNetwork{col(i)}, :, :);
        resultMatrixLogical = resultMatrix > 0;
        
        networkPairBlob = [];
        for j = 1:size(networkPair,1)
            for k = 1:size(networkPair,2)  
                networkPairBlob(j,k) = any(any(single(squeeze(networkPair(j,k,:,:))) .* resultMatrixLogical));
            end
        end
        
        %Place connections back into full matrix so they can be referenced
        %when pulling data from subjects
        fullMat = zeros(size(tStatsSignAll,1),size(tStatsSignAll,2));
        fullMat(parcelsIdxNetwork{row(i)}, parcelsIdxNetwork{col(i)}) = networkPairBlob;
        
        %Find connection idx in matrix
        connIdx = find(fullMat > 0);
        
        %Go through subject files and time-freq pull data for these connections
        %
        blobConnMatGroup1 = zeros(length(connIdx), size(resultMatrix,1), size(resultMatrix,2), size(group1Files,2));
        blobConnMatGroup2 = zeros(length(connIdx), size(resultMatrix,1), size(resultMatrix,2), size(group2Files,2));
        
        for n = 1:size(group1Files, 2)
            load(horzcat(group1Folder, group1Files{n}), '-mat');
            
            %This fcn will get the full list of idx of the stored data
            nonzeroIdxFull = reshape(bsxfun(@plus, double(subjectEdgeIdx), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]),[],1);
            %Then will re-create the original large matrix (n x n x freq x time)
            connMatFull = zeros(dims);
            connMatFull(nonzeroIdxFull) = nonzeroConn;
            connMatFullReshape = reshape(connMatFull, [size(connMatFull,1)*size(connMatFull,2), size(connMatFull,3), size(connMatFull,4)]);
            blobConnMatGroup1(:,:,:,n) = connMatFullReshape(connIdx,:,:);
            
            disp(['Loaded blob connections for Group1 subject ', num2str(n), '/', num2str(size(group1Files, 2))])
        end
        
        for n = 1:size(group2Files, 2)
            load(horzcat(group2Folder, group2Files{n}), '-mat');
            
            %This fcn will get the full list of idx of the stored data
            nonzeroIdxFull = reshape(bsxfun(@plus, double(subjectEdgeIdx), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]),[],1);
            %Then will re-create the original large matrix (n x n x freq x time)
            connMatFull = zeros(dims);
            connMatFull(nonzeroIdxFull) = nonzeroConn;
            connMatFullReshape = reshape(connMatFull,[size(connMatFull,1)*size(connMatFull,2),size(connMatFull,3),size(connMatFull,4)]);
            blobConnMatGroup2(:,:,:,n) = connMatFullReshape(connIdx,:,:);
            
            disp(['Loaded blob connections for Group2 subject ', num2str(n), '/', num2str(size(group2Files, 2))])
        end
        %}
        
        %For each connection, do t-test as usual to detect
        %within-connection blobs and find more significant connections
        pValuesReshape = []; tStatsReshape = [];
        numTests = 1000; permDist = zeros(numTests,length(connIdx));
        
        for n = 1:length(connIdx)
            %Pull data for connection
            conn1 = squeeze(blobConnMatGroup1(n,:,:,:));
            conn2 = squeeze(blobConnMatGroup2(n,:,:,:));
            
            %Baseline data
            meanConn1 = mean(conn1(:,baselineIdx,:),2);
            meanConn2 = mean(conn2(:,baselineIdx,:),2);
            
            conn1Baselined = conn1 - meanConn1;
            conn2Baselined = conn2 - meanConn2;
            %Reshape
            conn1Reshape = reshape(conn1Baselined, [size(conn1,1)*size(conn1,2),size(conn1,3)]);
            conn2Reshape = reshape(conn2Baselined, [size(conn2,1)*size(conn2,2),size(conn2,3)]);
            
            %Remove subs with zero connectivity
            conn1ReshapeNoZero = conn1Reshape(:,any(conn1Reshape));
            conn2ReshapeNoZero = conn2Reshape(:,any(conn2Reshape));
            
            [~, pValues, ~, stats] = ttest2(conn1ReshapeNoZero', conn2ReshapeNoZero');

            %Store results in this network pair matrix
            pValuesReshape(n, :, :) = reshape(pValues,     [size(conn1,1) size(conn1,2)]);
            tStatsReshape(n, :, :)  = reshape(stats.tstat, [size(conn1,1) size(conn1,2)]);   
            
            figure; imagesc(squeeze(pValuesReshape(n,:,:))<0.05)
            
            %%      Permutation testing on edges
            connAllReshape = [conn1ReshapeNoZero conn2ReshapeNoZero];

            for permNum = 1:numTests
                %Permute IDs
                permIdx = randperm(size(connAllReshape,2));

                %Regroup data
                permGroup1Conn = connAllReshape(:, permIdx(1:size(conn1ReshapeNoZero, 2)));
                permGroup2Conn = connAllReshape(:, permIdx(size(conn1ReshapeNoZero, 2)+1:end));

                %Perform test
                [~, pValues, ~, stats] = ttest2(permGroup1Conn', permGroup2Conn');
                
                pValuesPerm = reshape(pValues,     [size(conn1,1) size(conn1,2)]);
                tStatsPerm  = reshape(stats.tstat, [size(conn1,1) size(conn1,2)]); 
                
                %Find biggest blob
                tStatsPermMasked = tStatsPerm.*(pValuesPerm < 0.05);
                
                blobArrayTemp = [];
                [labeledMat, num] = bwlabel(tStatsPermMasked);
                if num > 0
                    for label = 1:num
                        tStatsLabel = labeledMat == label;
                        blobArrayTemp = [blobArrayTemp sum(sum(tStatsPermMasked.*tStatsLabel))];
                    end

                    permDist(permNum, n) = max(abs(blobArrayTemp));
                end
                if mod(permNum,50) == 0
                    disp(['Finished perm iteration ', num2str(permNum),' for connection ',num2str(n), '/', num2str(length(connIdx))])
                end
            end
          
        end
        
        %Find 95th percentile
        permDistSorted = sort(permDist,2,'descend');
        threshold = prctile(permDistSorted(:,1), 95);

        % Find significant blobs in true data vs perm
        tStatsReshapeMasked = tStatsReshape.*(pValuesReshape < 0.05);
        for conn = 1:size(tStatsReshapeMasked,1)
            connDataMasked = squeeze(tStatsReshapeMasked(conn,:,:));
            
            [labeledMat, num] = bwlabel(connDataMasked);
            for label = 1:num
                tStatsLabel = labeledMat == label;
                
                blobSum = sum(sum(connDataMasked.*tStatsLabel));
                
                if abs(blobSum) < threshold
                    labeledMat(labeledMat==label) = 0;
                end
            end
            if sum(sum(labeledMat)) == 0
                disp('Connection eliminated')
            else
                disp('Connection survived')
            end
        end
    end
end



%Goal is to first find which connections contribute to this difference,
%then find the ones with the biggest differences
1;
%

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

%}



function edgeFull = plotNCA_Matrix(resultMatrix, rowIdx, colIdx, divisionLineIdx, tStatsSignAll, networkNames, frequencies, latencies)
    
    %Find "middle" time-freq bin of blob to plot connections for
    %To account for non-circular blobs, will find the bin which is minimal
    %distance to all other bins (allowing for weird shaped blobs)
    [freqs, times] = find(resultMatrix > 0);
    D = sum(pdist2([freqs times], [freqs times],'euclidean'));
    middleIdx = find(D == min(D));
    xMid = freqs(middleIdx); yMid = times(middleIdx);
    
    %Plot surviving connections
    figure; 
    ax1 = subplot(1,2,1);
    imagesc(tStatsSignAll(:, :, xMid, yMid)); set(gcf, 'color', 'w')

    %Draw network division lines
    arrayfun(@(a)xline(a), divisionLineIdx(2:end-1))    
    arrayfun(@(a)yline(a), divisionLineIdx(2:end-1))

    %Plot diagonal
    xIdx = [1, size(tStatsSignAll,2)]; yIdx = [1, size(tStatsSignAll,1)];
    line(xIdx,yIdx,'color',[0 0 0])

    %Get indices for labeling plot (middle of network coordinates)
    labelCoords = round(mean([divisionLineIdx(1:end-1);divisionLineIdx(2:end)]));
    xticks(labelCoords); xticklabels(networkNames); xtickangle(45);
    yticks(labelCoords); yticklabels(networkNames); 
    set(gca,'TickLength',[0 0]);

    %Color shade connections based on direction of group difference
    colormap(ax1, [0 0 1 ; 1 1 1; 1 0 0])
    set(gca,'XAxisLocation','top','YAxisLocation','left');

    %Make matrix square (visualization purpose)
    x0 = 100; y0 = 100; width = 1000; height = 800; 
    set(gcf,'position',[x0, y0, width, height]); axis square
    xlabel('From Region'); ylabel('To Region')

    
    %Shading square
    parcelBounds1 = divisionLineIdx(rowIdx)+1:divisionLineIdx(rowIdx+1);
    parcelBounds2 = divisionLineIdx(colIdx)+1:divisionLineIdx(colIdx+1);
    
    group1SignCount(1) = length(find(tStatsSignAll(parcelBounds1, parcelBounds2, xMid, yMid) == 1));  %Higher conn in group1
    group2SignCount(2) = length(find(tStatsSignAll(parcelBounds1, parcelBounds2, xMid, yMid) == -1)); %Higher conn in group2    

    colorFactor = tanh(2*(group1SignCount - group2SignCount)/(group1SignCount + group2SignCount));
    
    boxColor = [0.5 + 0.5*colorFactor , 0, 0.5 - 0.5*colorFactor];
    patch([divisionLineIdx(colIdx), divisionLineIdx(colIdx), divisionLineIdx(colIdx+1), divisionLineIdx(colIdx+1)],... %x coordinates (box)
          [divisionLineIdx(rowIdx), divisionLineIdx(rowIdx+1), divisionLineIdx(rowIdx+1), divisionLineIdx(rowIdx)],... %y coordinates (box)
          boxColor, 'FaceAlpha', 0.4)    
      
    annotation('textbox', [0.45, 0.01, 0.1, 0.1], 'string',...
        horzcat('Extracted from Time: ',num2str(latencies(min(times))), ' to ', num2str(latencies(max(times))),'s, Freq: ',...
                                        num2str(frequencies(min(freqs))), ' to ', num2str(frequencies(max(freqs)))));
    
                                    
    %Plot time-frequency blob of results
    subplot(3,2,4);
    imagesc(latencies, 1:length(frequencies), resultMatrix > 0, [-1, 1]); 
    axis xy; hold on
    contour(latencies, 1:length(frequencies), resultMatrix > 0, 1, 'color', 'black', 'linewidth', 2.5) 
    YTickInterval = round(length(frequencies)/8);
    YTickIdx      = YTickInterval:YTickInterval:length(frequencies);
    YTickLabel    = round(round(frequencies(YTickIdx)*10)/10);
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickLabel);                                
                                    
    ylabel('Frequency (Hz)')
    set(get(gca,'YLabel'),'Rotation',90)  
    
    xlabel('Time (s)')
    
                            
    %Create Edge matrix for connectivity plot in BrainNet viewer
    edgeFull = zeros(size(tStatsSignAll,1),size(tStatsSignAll,2));
    edgeFull(divisionLineIdx(rowIdx)+1:divisionLineIdx(rowIdx+1),divisionLineIdx(colIdx)+1:divisionLineIdx(colIdx+1)) = tStatsSignAll(divisionLineIdx(rowIdx)+1:divisionLineIdx(rowIdx+1),divisionLineIdx(colIdx)+1:divisionLineIdx(colIdx+1), xMid, yMid);
    edgeFull(edgeFull == -1) = 0.5;

end

function [pValuesNetwork, tStatsNetwork] = runStatistics(conn1, conn2)
    %Run pairwise test on edges within network
    
    %Array Setup
    pValuesNetwork = single(ones(size(conn1,1), size(conn1,2), size(conn1,3)));
    tStatsNetwork  = single(zeros(size(conn1,1), size(conn1,2), size(conn1,3)));

    for edge = 1:size(conn1, 1)
        %Select data for edge
        edgeData1 = squeeze(conn1(edge, :, :, :));
        edgeData2 = squeeze(conn2(edge, :, :, :));

        edgeData1Reshape = reshape(edgeData1, [size(edgeData1,1)*size(edgeData1,2), size(edgeData1,3)]);
        edgeData2Reshape = reshape(edgeData2, [size(edgeData2,1)*size(edgeData2,2), size(edgeData2,3)]);

        %Leave only subjects with data
        edgeData1Nonzero = edgeData1Reshape(:, any(edgeData1Reshape));
        edgeData2Nonzero = edgeData2Reshape(:, any(edgeData2Reshape));

        %t-test
        if ~(isempty(edgeData1Nonzero) || isempty(edgeData2Nonzero))
            [~, pValues, ~, stats] = ttest2(edgeData1Nonzero', edgeData2Nonzero');

            %Store results in this network pair matrix
            pValuesNetwork(edge, :, :) = reshape(pValues,     [size(edgeData1,1) size(edgeData1,2)]);
            tStatsNetwork(edge, :, :)  = reshape(stats.tstat, [size(edgeData1,1) size(edgeData1,2)]);    
        end
    end

end


