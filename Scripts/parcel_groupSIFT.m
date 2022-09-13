% This function will perform groupSIFT-like analysis on a particular "seed"
% parcel
clear
%%                          Parameter setup
seedParcel = 221; %ACC (of salience) is nodes 221 and 223 of reordered list, 'Frontal.Mid.R'; %159 FrontoParietalR8 'Frontal.Sup.R'; %102 DefaultR13
subjThreshold = 50;
pixelPValue = 0.05;
permNum = 10000;
clusterLevelPvalue = 0.05;
GFWER = 1;

%%                      Selecting files to Load
%
%Select group 1 data
[group1Files, group1Folder] = uigetfile('*.conn', 'MultiSelect', 'on');
if ~any(group1Folder); disp('Cancelled.'); return; end

%Select group 2 data
[group2Files, group2Folder] = uigetfile('*.conn', 'MultiSelect', 'on');
if ~any(group1Folder); disp('Cancelled.'); return; end

tableGordon = table2cell(readtable('D:\NetworkConnectivity\Gordon2016_Map\Gordon286Nodes.node', 'filetype', 'text'));
%
%%                          
%First get network indices
load(horzcat(group1Folder, group1Files{1}), '-mat', 'parcelInfo', 'frequencies', 'latencies', 'dims');

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


%Load in group 1 data and pull out toIdx x fromIdx x t x f matrix
[connMatSegGroup1,finallySelectedEdgeIdx1, ~] = loadConnFiles(group1Files, group1Folder, seedParcel, dims, frequencies, latencies);

%Load in group 2 data and pull out toIdx x fromIdx x t x f matrix
[connMatSegGroup2,finallySelectedEdgeIdx2, dims] = loadConnFiles(group2Files, group2Folder, seedParcel, dims, frequencies, latencies);


%

%%              Calculate connections meet subject threshold

group1Subj = sum(squeeze(sum(sum(connMatSegGroup1,3),4))~=0,3);
group1Pct = group1Subj(1,:)/size(connMatSegGroup1,5);
group1Conns = find(100*group1Pct >= subjThreshold);

group2Subj = sum(squeeze(sum(sum(connMatSegGroup2,3),4))~=0,3);
group2Pct = group2Subj(1,:)/size(connMatSegGroup2,5);
group2Conns = find(100*group2Pct >= subjThreshold);

meetThreshold = intersect(group1Conns, group2Conns);


connGroup1Thresh = connMatSegGroup1(:,meetThreshold,:,:,:);
connGroup2Thresh = connMatSegGroup2(:,meetThreshold,:,:,:);

%%                      True t-statistics

for d = 1:2 %bidirectional, d=1 is TO seed, d=2 is FROM seed
    for i = 1:size(meetThreshold,2)
        connGroup1 = squeeze(connGroup1Thresh(d,i,:,:,:));
        connGroup2 = squeeze(connGroup2Thresh(d,i,:,:,:));
        
        if ndims(connGroup1) < 3
            connGroup1 = reshape(connGroup1, [size(connGroup1, 1), 1, size(connGroup1, 2)]);
            connGroup2 = reshape(connGroup2, [size(connGroup2, 1), 1, size(connGroup2, 2)]);
        end

        %Keep only subjects with nonzero values
        connGroup1Nonzero = connGroup1(:, :, any(squeeze(sum(sum(connGroup1, 2), 1)),2));
        connGroup2Nonzero = connGroup2(:, :, any(squeeze(sum(sum(connGroup2, 2), 1)),2));
        
        [edgeBlobMask, edgeTScore, edgePValue, edgeSurroMassOfCluster] = clusterLevelPermutationTest_oneside(connGroup1Nonzero, connGroup2Nonzero, 0, pixelPValue, permNum);
 
        clusterMask(d, i, :, :)    = edgeBlobMask;
        tStatistics(d, i, :, :)    = edgeTScore;
        pValues(d, i, :, :)    = edgePValue;
        surroMassOfCluster(:, d, i) = edgeSurroMassOfCluster(:);     
        disp(sprintf('edge %d complete', i))   
    end
end
%}
surroMassOfClusterReorg = [squeeze(surroMassOfCluster(:,1,:)) squeeze(surroMassOfCluster(:,2,:))];
surroMassOfClusterSorted = sort(surroMassOfClusterReorg,2,'descend');
if GFWER == 0
    criticalMassOfCluster = prctile(surroMassOfClusterSorted(:,1), 100-clusterLevelPvalue*100);
elseif GFWER == 1
    criticalMassOfCluster = prctile(surroMassOfClusterSorted(:,2), 100-clusterLevelPvalue*100);
end

%% Blob Correction (Keep only those > critical mass)
for d = 1:2
    for n = 1:size(meetThreshold,2)
        %Find surviving Blobs
        tmpBlobMask    = squeeze(clusterMask(d,n,:,:));
        tmpTStatistics = squeeze(tStatistics(d,n,:,:));
        [entryCount, blobId]  = hist(tmpBlobMask(:), unique(tmpBlobMask(:)));

        %Check each blob for significance
        for b = 2:length(blobId)
            blobMask = tmpBlobMask == blobId(b);
            blobSum = sum(sum(abs(tmpTStatistics .* blobMask)));

            if blobSum < criticalMassOfCluster
            %Remove blob if not significant
                tmpBlobMask = tmpBlobMask .* (~blobMask);
            else
                disp(['Surviving blob found in cell ' num2str(d),', ', num2str(n)]);
                
                %Find which parcels to plot
                if d == 1
                    toParcel = seedParcel;
                    fromParcel = meetThreshold(n);
                else
                    toParcel = meetThreshold(n);
                    fromParcel = seedParcel;
                end

                % set 2 dipoles
                sources(1,1).posxyz = str2num(parcelInfo{fromParcel,5});
                sources(1,1).momxyz = [0 0 0];
                sources(1,1).rv     = 0;
                sources(1,2).posxyz = str2num(parcelInfo{toParcel,5});
                sources(1,2).momxyz = [0 0 0];
                sources(1,2).rv     = 0;
                sources(1,3).posxyz = str2num(parcelInfo{toParcel,5}); %Duplicated because otherwise won't plot 2nd source
                sources(1,3).momxyz = [0 0 0];
                sources(1,3).rv     = 0;

                % 2 colors names officially supported by W3C specification for HTML
                colors{1,1}  = [1.00 0.66 0.76]; % source == red
                colors{2,1}  = [0.66 0.76 1.00]; % sink == blue

                % define colors used for dipoles
                plotColor = {colors{1,1} colors{2,1}};
                
                figure;
                coords = [0 90; 90 0; 0 0];
                for viewnum = 1:3
                    subplot(2,3,viewnum)
                    
                    %Plot dipoles
                    dipplot(sources, 'projimg', 'on', 'color', plotColor,  'dipolesize', 40,...
                                     'projlines', 'off', 'coordformat', 'MNI', 'spheres', 'on', 'gui', 'off');
                    set(findall(gca, '-property', 'linewidth'), 'linewidth', 3);       
                    view([coords(viewnum,:)]);
                    hold on 
                    
                    %Add lines that connects dipoles
                    posxyz(1,:) = sources(1,1).posxyz;
                    posxyz(2,:) = sources(1,2).posxyz;
                    arrow3d([posxyz(1,1) posxyz(2,1)], [posxyz(1,2) posxyz(2,2)], [posxyz(1,3) posxyz(2,3)], 0.7, 3, 7, [0.66 0.66 0.34]);
                end

                %Grab significant values/plot
                group1Values = squeeze(connGroup1Thresh(d,n,:,:,:));
                group2Values = squeeze(connGroup2Thresh(d,n,:,:,:));

                if size(tmpBlobMask,2) == 1 
                  %Line graph for resting state  
                    %Remove zeros for plotting
                    if ndims(group1Values) < 3
                        group1Values = reshape(group1Values, [size(group1Values, 1), 1, size(group1Values, 2)]);
                        group2Values = reshape(group2Values, [size(group2Values, 1), 1, size(group2Values, 2)]);
                    end

                    %Keep only subjects with nonzero values
                    group1ValuesNonzero = squeeze(group1Values(:, :, any(squeeze(sum(sum(group1Values, 2), 1)),2)));
                    group2ValuesNonzero = squeeze(group2Values(:, :, any(squeeze(sum(sum(group2Values, 2), 1)),2)));
                    
                    %Take mean of blob
                    meanGroup1 = mean(group1ValuesNonzero, 2);
                    meanGroup2 = mean(group2ValuesNonzero, 2);
                    
                    hold on
                    subplot(2,1,2)
                    C = {'k','b','r','g',};
                    legendLabels = {'Active', 'Sham'};
                
                    plot(meanGroup1, '.-', 'linewidth', 2,'color','r')
                    hold on
                    plot(meanGroup2, '.-', 'linewidth', 2,'color','b')
                
                    xlabel('Freq (Hz)'); ylabel('Connectivity');
                    title(['Active vs Sham: parcel ',num2str(fromParcel), '/', tableGordon{fromParcel,6}, ' (MNI: ', parcelInfo{fromParcel,5}, ')',...
                                             ' to ', num2str(toParcel), '/', tableGordon{toParcel,6}, ' (MNI: ', parcelInfo{toParcel,5}, ')'])
                    grid on
                    ylim([0, 1.2*max([meanGroup1;meanGroup2])])
                        
                    %Plotting transparent significance windows
                    blobLims = find(blobMask == 1);
                    blobBounds = [min(blobLims) max(blobLims)]; 
                    for b = 1:size(blobBounds,1)
                        xline(blobBounds(b,1),'r');
                        xline(blobBounds(b,2),'r');
                        %[blobMask(b,1) minY-0.01*minY blobBounds(b,2)-blobBounds(b,1) ((maxY+0.01*maxY)-(minY-0.01*minY))]
                        rectangle('Position',... % draw a semi-transparent rectangle
                          [blobBounds(b,1) 0 blobBounds(b,2)-blobBounds(b,1) 1.2*max([meanGroup1;meanGroup2])],...
                          'FaceColor',[1 0 0 0.1],'EdgeColor','none');
                    end
                    legend(legendLabels)
                    set(gcf, 'color', [1 1 1]);
                    
                    %For EXCEL export
                    group1ValuesMean = mean(squeeze(group1Values(blobMask,:,:)),1)';
                    group2ValuesMean = mean(squeeze(group2Values(blobMask,:,:)),1)';
                    
                    subjArray = [cellfun(@(x)x(1:3),group1Files','UniformOutput',false);cellfun(@(x)x(1:3),group2Files','UniformOutput',false)];
                    dataArray = [group1ValuesMean;group2ValuesMean];
                    1;
                else
                  %TIme-freq
                  
                end
                
            end
        end
        
        %Check if task or resting state data
        if size(tmpBlobMask,1) ~= 1
            adjustedMask(d, n, :, :) = tmpBlobMask;
        else %Rest
            adjustedMask(d, n, :) = tmpBlobMask;
        end
    end
end
    
1;




function [connMatSegGroup, finallySelectedEdgeIdx, dims] = loadConnFiles(groupFiles, groupFolder, seedParcel, dims, frequencies, latencies)    
    %Go through group files and extract connectivity for a network pair
    
    connMatSegGroup     = single(zeros(2, dims(2), length(frequencies), length(latencies), size(groupFiles,2)));
    connMatSegGroupTemp = single(zeros(2, dims(2), length(frequencies), length(latencies)));
    fprintf('[Group] Loaded file ')
    
    for n = 1:size(groupFiles, 2)
        load(horzcat(groupFolder, groupFiles{n}), '-mat');

        %This fcn will get the full list of idx of the stored data
        nonzeroIdxFull = reshape(bsxfun(@plus, double(subjectEdgeIdx), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]),[],1);
  
        %Then will re-create the original large matrix (n x n x freq x time)
        connMatFull = zeros(dims);
        connMatFull(nonzeroIdxFull) = nonzeroConn;
        
        toSeedMatrix = connMatFull(seedParcel, :, :, :);
        connMatSegGroup(1,:,:,:,n) = squeeze(toSeedMatrix);
        
        fromSeedMatrix = connMatFull(:, seedParcel, :, :);
        connMatSegGroup(2,:,:,:,n) = squeeze(fromSeedMatrix);
        
        
        %Can sanity check that new option above is correct
        %isequal(connMatSegGroupTemp, connMatFull(toIdx, fromIdx, :, :))
        
      %  disp(horzcat('[Group] Loaded file ', num2str(n), ' for network combo: ', num2str(network1), ' ', num2str(network2)))
        fprintf(['..',num2str(n)])
    end
    fprintf('\nFinished loading group files.\n')
end
