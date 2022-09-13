%
%%                      Selecting files to Load
%
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
numTests = 1000; %number of permutaions iterations, 200 to test but could be more 1000 for real analysis
countByNetwork = int8(zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), length(frequencies), length(latencies))); %This stores counts of each network for each time/freq bin
countByNetworkPerm = int8(zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), length(frequencies), length(latencies), numTests)); %Store of counts for perm tests

pValuesAll = uint8(ones(totalParcels, totalParcels, length(frequencies), length(latencies))); %This stores all p values for entire matrix (for plotting mainly)
%tStatsSignAll = int8(zeros(totalParcels,totalParcels,length(frequencies), length(latencies))); %This stores signs for all matrix

countByNetworkMaxStatistic = zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), numTests);

for network1 = 1:size(parcelsIdxNetwork, 2)
    for network2 = 1:size(parcelsIdxNetwork, 2)
        %%                      Load in Data
        toIdx   = parcelsIdxNetwork{network1};
        fromIdx = parcelsIdxNetwork{network2};
        
        %Load in group 1 data and pull out toIdx x fromIdx x t x f matrix
        [connMatSegGroup1,finallySelectedEdgeIdx1, ~] = loadConnFiles(group1Files, group1Folder, toIdx, fromIdx, frequencies, latencies, network1, network2);
        
        %Load in group 2 data and pull out toIdx x fromIdx x t x f matrix
        [connMatSegGroup2,finallySelectedEdgeIdx2, dims] = loadConnFiles(group2Files, group2Folder, toIdx, fromIdx, frequencies, latencies, network1, network2);

        %%                      Baseline data
        if dims(4) > 1
            %baseline = [-3.5 -2.5];
            baseline = [-1.6 -0.5];  %in seconds
            baselineIdx = find(latencies > baseline(1) & latencies < baseline(2));

            %Group 1 baseline
            baseline1   = mean(connMatSegGroup1(:, :, :, baselineIdx,:),4);
            connMatSegGroup1 = connMatSegGroup1 - baseline1;

            %Group 2 baseline
            baseline2   = mean(connMatSegGroup2(:, :, :, baselineIdx,:),4);
            connMatSegGroup2 = connMatSegGroup2 - baseline2;

            disp('Connectivity data baselined')
            clear baseline1 baseline2
        else
            baselineIdx = NaN;
        end
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
save(['/u/project/alenarto/gmicheli/CIDAR-GENETICS/allsubjs/reprocessed/15ics/' 'ADHD_Controls_test_tStatistics'],...
    'countByNetwork', 'tStatsSignAll', 'countByNetworkPerm', 'dims', ...
    'numTests', 'networkNames', 'frequencies', 'latencies', 'parcelInfo',...
    'parcelsIdxNetwork', 'divisionLineIdx', 'baselineIdx', ...
    'group1Folder', 'group1Files', 'group2Folder', 'group2Files', '-v7.3');
%}
1
%Save the above to a file for usage on diff fwer control p values...

%
%%          Now do FWER control and check out significant results
controlledNetworkResults = zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), length(frequencies), length(latencies));
countByNetworkMaxStatistic = zeros(size(parcelsIdxNetwork, 2), size(parcelsIdxNetwork, 2), size(countByNetworkPerm,5));

%FWER control in each network
for network1 = 1:size(parcelsIdxNetwork, 2)
    for network2 = 1:size(parcelsIdxNetwork, 2)
        
        countMatrix     = single(squeeze(countByNetwork(network1, network2, :, :)));
        countMatrixPerm = single(squeeze(countByNetworkPerm(network1, network2, :, :, :)));
        
        %Reshape for resting state if needed
        if ndims(countMatrixPerm) == 2 
            countMatrixPerm = reshape(countMatrixPerm,[size(countMatrixPerm,1),1,size(countMatrixPerm,2)]);
        end
        
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
    signifCorrNetwork = zeros(dims(1), dims(2));
    for i = 2:length(row)
        %%          Plotting significant networks & blobs
        resultMatrix = squeeze(controlledNetworkResults(row(i), col(i), :, :));
        edgeFull = plotNCA_Matrix(resultMatrix, row(i), col(i), divisionLineIdx, tStatsSignAll, networkNames, frequencies, latencies);
      
        %%               Load Subj Data for Network Pair
        %
        toIdx = parcelsIdxNetwork{row(i)};
        fromIdx = parcelsIdxNetwork{col(i)};
        
        %Load group 1 for network pair
        [connMatGroup1, finallySelectedEdgeIdx1] = loadConnFiles(group1Files, group1Folder, toIdx, fromIdx, frequencies, latencies, row(i), col(i));
        
        %Load group 2 for network pair
        [connMatGroup2, finallySelectedEdgeIdx2] = loadConnFiles(group2Files, group2Folder, toIdx, fromIdx, frequencies, latencies, row(i), col(i));
        
        if dims(4) > 1
            %Group 1 baseline
            baseline1   = mean(connMatGroup1(:, :, :, baselineIdx,:),4);
            connMatGroup1 = connMatGroup1 - baseline1;

            %Group 2 baseline
            baseline2   = mean(connMatGroup2(:, :, :, baselineIdx,:),4);
            connMatGroup2 = connMatGroup2 - baseline2;
        end
        %Select out overlapping edges (as we did before)
        overlapEdge = intersect(finallySelectedEdgeIdx1, finallySelectedEdgeIdx2);  
        overlapEdgeMat = zeros(dims(1), dims(2));
        overlapEdgeMat(overlapEdge) = 1;
        connMatGroup1 = connMatGroup1 .* overlapEdgeMat(toIdx, fromIdx);
        connMatGroup2 = connMatGroup2 .* overlapEdgeMat(toIdx, fromIdx);
        %This gives us network matrix for all valid connections

        %}
        %%        Choose which connections to run NCA regression on
        % Options:
        %   1. Run regression NCA on entire original network
        %   2. Run regression NCA on only significant blob connections
        %   - Both give a conn x time x freq x subj matrix
        
        regressConnMethod = 'Blob'; %'Original' or 'Blob'
        
        switch regressConnMethod
            case 'Original'
                %Get idx of group-overlapping edges
                selectedEdges = find(overlapEdgeMat(toIdx, fromIdx) == 1);
        
            case 'Blob'
                %Get significant connections from tStatsSignAll
                networkPair = tStatsSignAll(toIdx, fromIdx, :, :);
                
                %Get blob
                resultMatrixLogical = resultMatrix > 0;
                
                %Find connections that are significant at some time in blob
                networkPairBlob = [];
                for j = 1:size(networkPair,1)
                    for k = 1:size(networkPair,2)  
                        networkPairBlob(j,k) = any(any(single(squeeze(networkPair(j,k,:,:))) .* resultMatrixLogical));
                    end
                end
                
                %Find connection idx in matrix
                selectedEdges = find(networkPairBlob == 1);

                clear networkPair resultMatrixLogical
        end
        
        %Reshape to conn x time x freq x subj for easier NCA
        connMatGroup1Reshape = reshape(connMatGroup1, [size(connMatGroup1,1)*size(connMatGroup1,2), size(connMatGroup1,3), size(connMatGroup1,4), size(connMatGroup1,5)]);
        connMatGroup2Reshape = reshape(connMatGroup2, [size(connMatGroup2,1)*size(connMatGroup2,2), size(connMatGroup2,3), size(connMatGroup2,4), size(connMatGroup2,5)]);

        connMatGroup1Selected = connMatGroup1Reshape(selectedEdges,:,:,:);
        connMatGroup2Selected = connMatGroup2Reshape(selectedEdges,:,:,:);     
        clear connMatGroup1Reshape connMatGroup2Reshape
        %%      Choose which time/freq to run NCA regression on 
        % Options:
        %   1. Run on entire time-frequency matrix
        %       - Keeps edge x time x freq x subj matrix
        %   2. Run on mean connectivity in detected time-freq blob 
        %       - Gives new edge x avg conn x subj matrix
        
        regressTFMethod = 'Blob'; %'Full' or 'Blob'
        
        switch regressTFMethod
            case 'Full'
                continue
            case 'Blob'
                %Find avg time-freq connectivity in blob for all connection
                connMatGroup1Avged = zeros(size(connMatGroup1Selected,1),size(connMatGroup1Selected,4));
                connMatGroup2Avged = zeros(size(connMatGroup2Selected,1),size(connMatGroup2Selected,4));
                
                for edgeNum = 1:length(selectedEdges)
                    %Group 1
                    edgeConn = squeeze(connMatGroup1Selected(edgeNum,:,:,:));
                    %Reshape for resting state if needed
                    if ndims(edgeConn) == 2 
                        edgeConn = reshape(edgeConn,[size(edgeConn,1),1,size(edgeConn,2)]);
                    end
                    edgeConnMasked = reshape(edgeConn, [size(edgeConn,1)*size(edgeConn,2), size(edgeConn,3)]);
                    connMatGroup1Avged(edgeNum,:) = mean(edgeConnMasked(find(resultMatrix~=0),:));
                    
                    %Group 2
                    edgeConn = squeeze(connMatGroup2Selected(edgeNum,:,:,:));
                    %Reshape for resting state if needed
                    if ndims(edgeConn) == 2 
                        edgeConn = reshape(edgeConn,[size(edgeConn,1),1,size(edgeConn,2)]);
                    end
                    edgeConnMasked = reshape(edgeConn, [size(edgeConn,1)*size(edgeConn,2), size(edgeConn,3)]);
                    connMatGroup2Avged(edgeNum,:) = mean(edgeConnMasked(find(resultMatrix~=0),:));
                end
                clear edgeConnMasked edgeConn
        end
        
        %%                  Create spreadsheet to export
        [~, ~, currentTable] = xlsread('D:/NetworkConnectivity/Gordon2016_Map/Parcels/Parcels_withNoneRemovedReOrdered.xlsx');
        
        %Find indices of significant regions
        [rowRegions, colRegions] = find(edgeFull~=0);
        uniqueRegions = unique([rowRegions;colRegions]);
        
        %Find corresponding idx in originalTable, so that .nii regions can
        %be accessed properly in BrainNet
        originalRegions = vertcat(parcelInfo{uniqueRegions,1})';

        % Build array of connectivity headers
        headerArray = {};
        for n = 1:length(rowRegions)
            headerArray{n} = [currentTable{colRegions(n),2},'_to_',currentTable{rowRegions(n),2}];
        end

        %%                  Pull in spreadsheet of data
        dataTable = readtable('D:/NetworkConnectivity/FlankerData/TicStudy_PerformanceMeasures.csv');
        
        %Code assumes connecticity data above is "stackable" in a way that
        %matches the spreadsheet

        %Will create connDataAll variable, which is a elements x subj
        %array, for each regression on each element
        switch regressTFMethod
            case 'Full'
                connDataAll = cat(4,connMatGroup1Selected,connMatGroup2Selected);
                connDataAll = reshape(connDataAll, [size(connDataAll,1)*size(connDataAll,2)*size(connDataAll,3), size(connDataAll,4)]);
            case 'Blob'
                connDataAll = cat(2,connMatGroup1Avged,connMatGroup2Avged);  
        end
        connDataAll(connDataAll==0) = NaN;
        
        %Run linear models on each element
        corrMethod = 'Corr'; %Corr, PartialCorr, or LinearModel
        
        pValueRegress = zeros(size(connDataAll,1),1);
        rhoRegress = zeros(size(connDataAll,1),1);
        for edge = 1:size(connDataAll,1)
            switch corrMethod
                case 'Corr'
                    [rhoRegress(edge), pValueRegress(edge)] = corr(connDataAll(edge, :)', dataTable.acc_valid_congruent_pct, 'rows', 'pairwise');
                case 'PartialCorr'
                    [rhoRegress(edge), pValueRegress(edge)] = partialcorr(connDataAll(edge, :)', dataTable.acc_valid_congruent_pct, dataTable.Age, 'rows', 'pairwise');
                case 'LinearModel'
                    dataTable.conn = connDataAll(edge, :)';
                    lm = fitlm(dataTable, 'conn ~ acc_valid_congruent_pct'); %Define model
                    pValueRegress(edge) = lm.Coefficients{'acc_valid_congruent_pct','pValue'};
            end
            if mod(edge,100) == 0
               disp(['Finished test regression on element ', num2str(edge), ' of ', num2str(size(connDataAll,1))])
            end
        end
        
        %Run permutation linear models
        numTests = 100;
        pValueRegressPerm = zeros(size(connDataAllPerm,1), numTests);
        for t = 1:numTests
            %Permute data
            permIdx        = randperm(size(connDataAll,2));
            connDataAllPerm = connDataAll(:, permIdx);

            %Run linear models on each element

            for edge = 1:size(connDataAllPerm,1)
                switch corrMethod
                    case 'Corr'
                        [rho, pValueRegress(edge)] = corr(connDataAllPerm(edge, :)', dataTable.acc_valid_congruent_pct, 'rows', 'pairwise');
                    case 'PartialCorr'
                        [rho, pValueRegress(edge)] = partialcorr(connDataAllPerm(edge, :)', dataTable.acc_valid_congruent_pct, dataTable.Age, 'rows', 'pairwise');
                    case 'LinearModel'
                        dataTable.conn = connDataAllPerm(edge, :)';
                        lm = fitlm(dataTable, 'conn ~ acc_valid_congruent_pct'); %Define model
                        pValueRegressPerm(edge, t) = lm.Coefficients{'acc_valid_congruent_pct','pValue'};
                end
            end
            
            if mod(t,10) == 0
               disp(['Finished test regression on iteration ', num2str(t), '/', num2str(numTests)])
            end
        end
        
        %Find & apply perm thresholds
        threshPValue = 0.05;
        switch regressTFMethod
            case 'Full'
                
            case 'Blob'
                pValueRegressPerm = sum(pValueRegressPerm < threshPValue);
                permThresh = prctile(pValueRegressPerm, 95);
                
                trueSignif = sum(pValueRegress < threshPValue);
                if trueSignif > permThresh
                    disp(['Significant Model found: Conn (To: ', networkNames{row(i)}, ', From: ', networkNames{col(i)}, ') ~ ', lm.PredictorNames{:}]);
                    
                    %Create indiv network matrix of corr coeffs
                    %Can either plot this by itself, or together with other
                    %significant network pairs in whole brain view
                    networkMatRegress = zeros(length(toIdx), length(fromIdx));
                    networkMatRegress(selectedEdges(pValueRegress < threshPValue)) = rhoRegress(pValueRegress < threshPValue);
                    
                    %Place into whole network
                    signifCorrNetwork(toIdx, fromIdx) = networkMatRegress;
                end
        end
    end
    
    %Plot NCA summary for regression with chosen variable
   % if any(any(signifCorrNetwork))
   %     plotNCA_Matrix_Regress(1, rowIdx, colIdx, divisionLineIdx, signifCorrNetwork, networkNames, 1, 1)
        
  %  end
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

function [connMatSegGroup, finallySelectedEdgeIdx, dims] = loadConnFiles(groupFiles, groupFolder, toIdx, fromIdx, frequencies, latencies, network1, network2)    
    %Go through group files and extract connectivity for a network pair
    
    connMatSegGroup = single(zeros(length(toIdx), length(fromIdx), length(frequencies), length(latencies), size(groupFiles,2)));
    connMatSegGroupTemp = single(zeros(length(toIdx), length(fromIdx), length(frequencies), length(latencies)));
    fprintf('[Group] Loaded file ')
    for n = 1:size(groupFiles, 2)
        tic
        load(horzcat(groupFolder, groupFiles{n}), '-mat');

        %This fcn will get the full list of idx of the stored data
        nonzeroIdxFull = reshape(bsxfun(@plus, double(subjectEdgeIdx), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]),[],1);
        
  %%%%%%%%%%%%%%% 8-2-2022 NEW - Helps files load faster
        %Get connection idx for selected network
        connMatTemp = zeros(dims(1),dims(2));
        connMatTemp(subjectEdgeIdx) = subjectEdgeIdx;
        connMatNetwork = connMatTemp(toIdx, fromIdx);
        connInNetwork = connMatNetwork(connMatNetwork>0);
        
        %Find idx of where time-freq conn values are in data file
        nonzeroIdxNetwork = reshape(bsxfun(@plus, double(connInNetwork), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]),[],1);
        [sharedvals,idx] = intersect(nonzeroIdxFull,nonzeroIdxNetwork,'stable');
        
        smallConnInNetwork = find(connMatNetwork>0);
        smallNonzeroIdxNetwork = reshape(bsxfun(@plus, double(smallConnInNetwork), size(connMatNetwork,1)*size(connMatNetwork,2)*[0:dims(3)*dims(4)-1]),[],1);
        connMatSegGroupTemp(smallNonzeroIdxNetwork) = nonzeroConn(idx);
        
        connMatSegGroup(:,:,:,:,n) = connMatSegGroupTemp;
  %%%%%%%%%%%%%%%      
  
        %Then will re-create the original large matrix (n x n x freq x time)
        %connMatFull = zeros(dims);
        %connMatFull(nonzeroIdxFull) = nonzeroConn;
        %connMatSegGroup(:,:,:,:,n) = single(connMatFull(toIdx, fromIdx, :, :));
        
        %Can sanity check that new option above is correct
        %isequal(connMatSegGroupTemp, connMatFull(toIdx, fromIdx, :, :))
        
      %  disp(horzcat('[Group] Loaded file ', num2str(n), ' for network combo: ', num2str(network1), ' ', num2str(network2)))
        fprintf(['..',num2str(n)])
        toc
    end
    fprintf('/nFinished loading group files./n')
end


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
    %imagesc(tStatsSignAll(:, :, xMid(1), yMid(1))); set(gcf, 'color', 'w')
    totalSignAll = sign(sum(sum(tStatsSignAll(:, :, freqs, times),3),4));
    imagesc(totalSignAll); set(gcf, 'color', 'w')
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
    if length(latencies)>1
        contour(latencies, 1:length(frequencies), resultMatrix > 0, 1, 'color', 'black', 'linewidth', 2.5) 
    end
    YTickInterval = round(length(frequencies)/8);
    YTickIdx      = YTickInterval:YTickInterval:length(frequencies);
    YTickLabel    = round(round(frequencies(YTickIdx)*10)/10);
    set(gca, 'YTick', YTickIdx, 'YTickLabel', YTickLabel);                                
                                    
    ylabel('Frequency (Hz)')
    set(get(gca,'YLabel'),'Rotation',90)  
    
    xlabel('Time (s)')
    
                            
    %Create Edge matrix for connectivity plot in BrainNet viewer
    edgeFull = zeros(size(tStatsSignAll,1),size(tStatsSignAll,2));
   % edgeFull(divisionLineIdx(rowIdx)+1:divisionLineIdx(rowIdx+1),divisionLineIdx(colIdx)+1:divisionLineIdx(colIdx+1)) = tStatsSignAll(divisionLineIdx(rowIdx)+1:divisionLineIdx(rowIdx+1),divisionLineIdx(colIdx)+1:divisionLineIdx(colIdx+1), xMid, yMid);
    edgeFull(divisionLineIdx(rowIdx)+1:divisionLineIdx(rowIdx+1),divisionLineIdx(colIdx)+1:divisionLineIdx(colIdx+1)) = totalSignAll(divisionLineIdx(rowIdx)+1:divisionLineIdx(rowIdx+1),divisionLineIdx(colIdx)+1:divisionLineIdx(colIdx+1));
    edgeFull(edgeFull == -1) = 0.5;

end


function [pValuesNetwork, tStatsNetwork] = runStatistics(conn1, conn2)
    %Run pairwise test on edges within network
    
    %Array Setup
    pValuesNetwork = single(ones(size(conn1,1), size(conn1,2), size(conn1,3)));
    tStatsNetwork  = single(zeros(size(conn1,1), size(conn1,2), size(conn1,3)));

    %Reshape for resting state if needed
    if ndims(pValuesNetwork) == 2 
        pValuesNetwork = reshape(pValuesNetwork,[size(pValuesNetwork,1),1,size(pValuesNetwork,2)]);
        tStatsNetwork = reshape(tStatsNetwork,[size(tStatsNetwork,1),1,size(tStatsNetwork,2)]);
    end
    
    for edge = 1:size(conn1, 1)
        %Select data for edge
        edgeData1 = squeeze(conn1(edge, :, :, :));
        edgeData2 = squeeze(conn2(edge, :, :, :));
        %Reshape for resting state if needed
        if ndims(edgeData1) == 2 
            edgeData1 = reshape(edgeData1,[size(edgeData1,1),1,size(edgeData1,2)]);
            edgeData2 = reshape(edgeData2,[size(edgeData2,1),1,size(edgeData2,2)]);
        end
        edgeData1Reshape = reshape(edgeData1, [size(edgeData1,1)*size(edgeData1,2), size(edgeData1,3)]);
        edgeData2Reshape = reshape(edgeData2, [size(edgeData2,1)*size(edgeData2,2), size(edgeData2,3)]);

        %Leave only subjects with data
        edgeData1Nonzero = edgeData1Reshape(:, any(edgeData1Reshape));
        edgeData2Nonzero = edgeData2Reshape(:, any(edgeData2Reshape));

        %t-test
        if ~(isempty(edgeData1Nonzero) || isempty(edgeData2Nonzero))
            try
            [~, pValues, ~, stats] = ttest2(edgeData1Nonzero', edgeData2Nonzero');
                        %Store results in this network pair matrix
            pValuesNetwork(edge, :, :) = reshape(pValues,     [size(edgeData1,1) size(edgeData1,2)]);
            tStatsNetwork(edge, :, :)  = reshape(stats.tstat, [size(edgeData1,1) size(edgeData1,2)]);    
            catch
            end

        end
    end

end


