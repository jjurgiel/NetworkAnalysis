clear
userInputFilename = 'TNSNetwork50p';
userInputPercentage = 0.5;

%% Loading in Data
[allFiles, workingFolder] = uigetfile('*.set', 'MultiSelect', 'on');
if ~any(workingFolder); disp('Cancelled.'); return; end

%Switch directories
cd(workingFolder)

% Load all .set data located under the working folder
ALLEEG = []; numICs = []; subIDList = {};
for n = 1:length(allFiles)
    EEG = pop_loadset('filename', allFiles{n},'filepath',workingFolder, 'loadmode', 'info');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    numICs(n) = size(EEG.icaweights,1);
    subIDList{n,1} = ALLEEG(n).subject;
end
fileNameList = {ALLEEG.filename}';

% Compute # of bidirection IC connections in each subject, minus self connections
numIcCombination = sum(numICs.^2);


%% Get unique dipole locations & connectivity data for each subject

fprintf('\nExtracting dipole locations & information for each subject...\n')

dimensionLabels= ALLEEG(1).CAT.Conn.dims;
dimensionLabels{1,5} = 'subjects';
timeFreqSize = size(squeeze(ALLEEG(1).CAT.Conn.RPDC(  1,1,:,:)));
latencies = ALLEEG(1).CAT.Conn.erWinCenterTimes;
frequencies = ALLEEG(1).CAT.Conn.freqs;

% dipolePairAndMeasureObj will store pairwise connectivity info for all subj
dipolePairAndMeasureObj                   = pr.dipolePairAndMeasure;
dipolePairAndMeasureObj.linearizedMeasure = zeros([numIcCombination timeFreqSize(1)*timeFreqSize(2)]);



dipoleCounter = 0;
counter = 1;
%Go through all selected EEG files
for n = 1:length(ALLEEG)
    numberOfSessiondipoles = length(ALLEEG(n).dipfit.model);
    dipoleLocation = [];
    dipoleResidualVariance = [];
    for m = 1:numberOfSessiondipoles
        dipoleLocation(m,:) = selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz);
        dipoleResidualVariance(m) = ALLEEG(n).dipfit.model(m).rv;
    end
    dipoleId = (dipoleCounter+1):(dipoleCounter+length(dipoleLocation));
    dipoleCounter = dipoleCounter + length(dipoleLocation);
    
    %Go through pairs of dipoles and extract connectivity for each edge
    for m = 1:length(dipoleLocation) % from location: confirm it with ALLEEG(1,1).CAT.Conn: dims: {'var_to'  'var_from'  'freq'  'time'}
        for k = 1:length(dipoleLocation) % to location
            dipolePairAndMeasureObj.from.location = cat(1, dipolePairAndMeasureObj.from.location, selectDipWithLargerMoment(ALLEEG(n).dipfit.model(m).posxyz, ALLEEG(n).dipfit.model(m).momxyz));
            dipolePairAndMeasureObj.to.location = cat(1, dipolePairAndMeasureObj.to.location, selectDipWithLargerMoment(ALLEEG(n).dipfit.model(k).posxyz, ALLEEG(n).dipfit.model(k).momxyz));
            
            if m == k
                dipolePairAndMeasureObj.linearizedMeasure(counter,:) = zeros(timeFreqSize(1)*timeFreqSize(2),1);
            else
                %Again, ALLEEG(1,1).CAT.Conn: dims: {'var_to'  'var_from'  'freq'  'time'}
                tmpTimeFreqMatrix = squeeze(ALLEEG(n).CAT.Conn.RPDC(  k,m,:,:));
                dipolePairAndMeasureObj.linearizedMeasure(counter,:) = tmpTimeFreqMatrix(:);
            end
            dipolePairAndMeasureObj.sessionId(counter,1) = n;
            dipolePairAndMeasureObj.fromDipoleId(counter,1) = dipoleId(m);
            dipolePairAndMeasureObj.toDipoleId(counter,1)   = dipoleId(k);
            counter = counter + 1;
        end
    end
end

% Find unique dipoles
uniqueDipole = pr.dipole;
[~,fromId,fromIdReverse] = unique(dipolePairAndMeasureObj.fromDipoleId, 'stable');

% Extract coordinates of the unique dipoles
uniqueDipole.location = dipolePairAndMeasureObj.from.location(fromId,:);


%% Define anatomical parcels
%Current design is set for Gordon mapping. If you want to use Yeo or
%something else, add Case below and confirm that coordinates of .nii files
%match those in EEGLAB
parcelMap = 'Gordon';
switch parcelMap
    case 'Yeo'
        niiInfo = niftiinfo( 'D:\NetworkConnectivity\1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels\MNI152\Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz');
        niiMatrix = niftiread('D:\NetworkConnectivity\1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels\MNI152\Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz');
        [~,~,networkTxt] = xlsread('D:\NetworkConnectivity\1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels\MNI152\Yeo2011_7networks_N1000.split_components.glossary.JJ.xlsx');
        xFlip = -1; %X coordinate (left/right) is backwards in Yeo
    case 'Gordon'
        niiInfo = niftiinfo( 'D:\NetworkConnectivity\Gordon2016_Map\Parcels\Parcels_MNI_222.nii');
        niiMatrix = niftiread('D:\NetworkConnectivity\Gordon2016_Map\Parcels\Parcels_MNI_222.nii');
      %  [~,~,networkTxt] = xlsread('D:\NetworkConnectivity\Gordon2016_Map\Parcels\Parcels.xlsx');
      %  [~,~,networkTxt] = xlsread('D:\NetworkConnectivity\Gordon2016_Map\Parcels\Parcels_withNoneRemoved.xlsx');
        [~,~,networkTxt] = xlsread('D:\NetworkConnectivity\Gordon2016_Map\Parcels\Parcels_withNoneRemovedReOrdered.xlsx');
        xFlip = 1;
end

%groupSIFT (89 x 106 x 89)
%        minGridLocation = [-88 -121  -77]; maxGridLocation = [89   90   99];
% Actual min for 2mm
%x -88 to 88, y  -121 to 89, z -77 to 99
%Parcel (YEO 7) (91 x 109 x 91) +90 -126 -72 [T matrix -2,2,2]
%Gordon (91 x 109 x 91) -90. -126, -72 (t matrix 2,2,2), same starting as YEO
%Actual brain coordinates for x/y/z min/max
%x/y/z = headGrid.zCube(headGrid.insideBrainCube);
%groupSIFT = x -70 to 72, y -107 to 75, z -61 to 85
%Parcel = x -70 to 70, y -107 to 71, 


fprintf('NOTE: .nii file Transformation matrix starting coordinates are: X = %d, Y = %d, Z = %d\n', ...
    niiInfo.Transform.T(4,1), niiInfo.Transform.T(4,2), niiInfo.Transform.T(4,3));
disp('EEGLAB assumes left hem X is negative, back of head Y negative, and bottom of head Z negative')
disp('This means none of the above should be positive. If positive, adjust code.')

if niiInfo.Transform.T(4,1) > 0
   niiInfo.Transform.T(4,1) = xFlip*niiInfo.Transform.T(4,1);
   niiInfo.Transform.T(1,1) = xFlip*niiInfo.Transform.T(1,1);
end

%% Set up some headgrid info & calculate dipole probabilities

%headgridSpacing set to 2mm to match .nii fMRI files
headGridSpacing = 2; 

tMat = niiInfo.Transform.T;
parcelHeadGrid = pr.headGrid(headGridSpacing,[tMat(4,1) tMat(4,2) tMat(4,3)],...
                 [tMat(4,1) + tMat(1,1)*(niiInfo.ImageSize(1)-1), ...
                  tMat(4,2) + tMat(2,2)*(niiInfo.ImageSize(2)-1), ...
                  tMat(4,3) + tMat(3,3)*(niiInfo.ImageSize(3)-1)]);
if any((size(parcelHeadGrid.insideBrainCube) ==  size(parcelHeadGrid.xCube))==0)
    error('insideBrainCube does not match dimensions of xCube, delete precomputed headGrid from +pr/precomputed')
end

% Create Gaussian Weight Matrix
standardDeviationOfEstimatedDipoleLocation = 20/2.355; % FWHM = 20, this calculates sigma in Gaussian equation
projectionParameter = pr.projectionParameter(standardDeviationOfEstimatedDipoleLocation);
projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian = 3;
[~,~,gaussianWeightMatrixParcel] = pr.meanProjection.getProjectionMatrix(uniqueDipole, parcelHeadGrid, projectionParameter);


parcelOfInterest          = [];
parcelNums                = vertcat(networkTxt{:,1});%length(unique(niiMatrix(:)))-1;
dipoleProbabilityInParcel = zeros(uniqueDipole.numberOfDipoles, length(parcelNums));

%Calculate parcel densities for each dipole to see where they project to
for i = 1:length(parcelNums)
    parcelOfInterest(i).label = networkTxt{i,1};
    parcelOfInterest(i).network = regexprep(networkTxt{i,2}, ' ', '_');
    parcelOfInterest(i).hemisphere = regexprep(networkTxt{i,3}, ' ', '_');
    parcelOfInterest(i).membershipProbabilityCube = permute(double((niiMatrix == parcelOfInterest(i).label)>0),[2,1,3]); 
    parcelOfInterest(i).probabilityThreshold = 0.001;
    %parcelOfInterest(i).membershipProbabilityCube = parcelOfInterest(i).membershipCube;
    
    %X COORD OF LEFT HEM PARCELS ARE ACTUALLY RIGHT HEM...
   dipoleProbabilityInParcel(:,i) = gaussianWeightMatrixParcel * parcelOfInterest(i).membershipProbabilityCube(parcelHeadGrid.insideBrainCube);
   
   %%%% SANITY CHECK PLOTTING FOR DIPOLE & PARCEL LOCATIONS %%%% 
%     parcelIdx = find(parcelOfInterest(i).membershipCube == 1);
%     xCubeTemp = xFlip*parcelHeadGrid.xCube(parcelIdx);
%     yCubeTemp = parcelHeadGrid.yCube(parcelIdx);
%     zCubeTemp = parcelHeadGrid.zCube(parcelIdx);
%     figure; hold on;
%     scatter3(xCubeTemp(:), yCubeTemp(:), zCubeTemp(:), 30, 'b')
%     for n = 1:length(dipoleProbabilityInParcel(:,i))
%         if dipoleProbabilityInParcel(n,i) ~= 0
%             hold on
%             scatter3(dipoleObj.location(n,1), dipoleObj.location(n,2),dipoleObj.location(n,3),30,'g')
%         else
%             scatter3(dipoleObj.location(n,1), dipoleObj.location(n,2),dipoleObj.location(n,3),30,'r')
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%
    fprintf('.')
end
fprintf('\nFinished estimating parcel probabilities.\n')


%%  Obtain subject % contributions to all parcels and find preselected ROI idx
icSubjIdList = [];
for subjId = 1:length(ALLEEG)
    tmpNumICs = ALLEEG(1,subjId).CAT.nbchan;
    icSubjIdList = [icSubjIdList; repmat(subjId, [tmpNumICs 1])];
end
subjDipoleDensityTable = zeros(length(ALLEEG), size(dipoleProbabilityInParcel,2));
for subjId = 1:length(ALLEEG)
    tmpSubjIdx = find(icSubjIdList == subjId);
    subjDipoleDensityTable(subjId,:) = sum(dipoleProbabilityInParcel(tmpSubjIdx,:),1);
end
clear ALLEEG %No longer needed
roiNonzeroSubjDetectionVector = sum(subjDipoleDensityTable ~= 0);

%%% Plot how many parcels meet subject threshold criteria %%%
subjCounts = histcounts(roiNonzeroSubjDetectionVector,20);
subjCumSum = [size(roiNonzeroSubjDetectionVector,2) size(roiNonzeroSubjDetectionVector,2)-cumsum(subjCounts)];
figure; bar(0:5:100,subjCumSum); xlabel('% subjects'); ylabel('Num parcels meeting % threshold');

%Select only ROIs which meet suject threshold criteria
preselectedRoiIdx = find(roiNonzeroSubjDetectionVector>= userInputPercentage*size(subjDipoleDensityTable,1)); 
fprintf('\nSubject percentages across parcels plotted.\n')

% Convert the preselected ROI indices to edge (dipolePairdensity) index.
% Max == N_parcel^2 - N_parcel (these are self referencing edges) 
roiLabelsReduced = vertcat(parcelOfInterest.label);
clear parcelOfInterest %No longer needed
preselectedRoiIdxRepeatedMatrix = repmat((preselectedRoiIdx'-1)*length(roiLabelsReduced), [1 length(preselectedRoiIdx')]);
preselectedEdgeIdx = bsxfun(@plus, preselectedRoiIdxRepeatedMatrix, preselectedRoiIdx)';
preselectedEdgeIdx(logical(eye(size(preselectedEdgeIdx)))) = 0; % Edge's self connection is excluded.
preselectedEdgeIdx = nonzeros(preselectedEdgeIdx(:)); %This is 1d array of matrix idx of connections


% Calculate dipole densities in source regions and destination regions. Reuse indices from 'from' dipole for 'to' dipole
fromDipoleProbabilityInRegion = zeros(dipolePairAndMeasureObj.from.numberOfDipoles, length(roiLabelsReduced));
toDipoleProbabilityInRegion   = zeros(dipolePairAndMeasureObj.to.numberOfDipoles,   length(roiLabelsReduced));

[~,~,toIdReverse] = unique(dipolePairAndMeasureObj.toDipoleId, 'stable');
[~,fromId,fromIdReverse] = unique(dipolePairAndMeasureObj.fromDipoleId, 'stable');
for i = 1:length(roiLabelsReduced)
    fromDipoleProbabilityInRegion(:,i) = dipoleProbabilityInParcel(fromIdReverse, i);
    toDipoleProbabilityInRegion(:,  i) = dipoleProbabilityInParcel(toIdReverse,   i);
end

% Prepare dipole pair density by multiplying from-density and to-density. The first dimension is the sum of each subject's IC*IC, the second is ROI^2
% Note that this maps the combination of ICs to ROI*ROI.
dipolePairDensityFromIcSquareToRoiSquare = zeros(dipolePairAndMeasureObj.from.numberOfDipoles, length(roiLabelsReduced) * length(roiLabelsReduced), 'single');
for i = 1:size(dipolePairDensityFromIcSquareToRoiSquare,1)
    tmpDipoleProbability = (fromDipoleProbabilityInRegion(i,:)' *  toDipoleProbabilityInRegion(i,:))';
    tmpDipoleProbability(logical(eye(size(tmpDipoleProbability)))) = 0; % Edge's self connection is excluded.
    dipolePairDensityFromIcSquareToRoiSquare(i, :) = vec(tmpDipoleProbability)';
end

% Compute dipole pair density for only preselected pairs (actively zero-outs the non-selected pairs; also saves time.)
numberOfsessions = max(dipolePairAndMeasureObj.sessionId); % Session means subjects.
dipolePairDensity = zeros(size(dipolePairDensityFromIcSquareToRoiSquare, 2), numberOfsessions);
for i = 1:length(preselectedEdgeIdx)
    currentEdgeIdx = preselectedEdgeIdx(i);
    currentEdgeIcSquare = dipolePairDensityFromIcSquareToRoiSquare(:,currentEdgeIdx);
    sessionDensity = zeros(1,numberOfsessions);
    for j=1:numberOfsessions
        sessionDensity(j) = sum(currentEdgeIcSquare(dipolePairAndMeasureObj.sessionId == j));
    end

    % Store dipole pair density in ROI^2 x subj.
    dipolePairDensity(currentEdgeIdx,:) = sessionDensity;
end

dipolePairDensity = reshape(dipolePairDensity, [length(roiLabelsReduced) length(roiLabelsReduced) numberOfsessions]);

%Find EDGES meeting subj threshold criteria (instead of later like
%groupsift)
%subjInEdges = sum(logical(dipolePairDensity),3);
%numSubjThreshold = userInputPercentage*size(subjDipoleDensityTable,1);
%groupEdgesMeetingThreshold = subjInEdges >= numSubjThreshold;

dipoleProbabilityInRegion = dipoleProbabilityInParcel;
linearizedMeasure = dipolePairAndMeasureObj.linearizedMeasure;
parcelInfo = networkTxt;
save([workingFolder userInputFilename '_dipolePairDensity'],...
    'dipolePairDensity', 'dipoleProbabilityInParcel', 'parcelInfo',...
    'fileNameList', 'dipolePairDensityFromIcSquareToRoiSquare', ...
    'preselectedRoiIdx', '-v7.3');

% prepare the map
numOfRegionsOfInterest = size(dipolePairDensity,1);
toRegionNumber   = repmat((1:numOfRegionsOfInterest)', [1, numOfRegionsOfInterest]); % first deimension (rows) contains toRegion Numbers.
toRegionNumber   = toRegionNumber(:);
fromRegionNumber = repmat( 1:numOfRegionsOfInterest,   [numOfRegionsOfInterest, 1]); % second deimension (columns) contains fromRegion Numbers.
fromRegionNumber = fromRegionNumber(:);

%allConnectivityStack = single(zeros(size(dipolePairDensity,1), size(dipolePairDensity,2), length(frequencies), length(latencies), length(allFiles)));
%allConnectivityStack = single(zeros(size(dipolePairDensity,1), size(dipolePairDensity,2),                   1, length(latencies), length(allFiles))); %No freq dimension
%allConnectivityStack = single(zeros(size(dipolePairDensity,1), size(dipolePairDensity,2),                   1,                 1, length(allFiles))); %No time/freq dimension

%First will get which connections meet subj threshold. This could be done
%cleaner but this is simplest way
countMatrix = zeros(max(fromRegionNumber), max(fromRegionNumber), length(unique(dipolePairAndMeasureObj.sessionId)'), 'single');
linearizedMeasureTemp = dipolePairAndMeasureObj.linearizedMeasure~=0;
for sessionIdx = unique(dipolePairAndMeasureObj.sessionId)'

    %Will use matrix multiplication instead of old way
    currentSessionIdx = find(dipolePairAndMeasureObj.sessionId==sessionIdx); 
    dipolePairDensitiesInIcSquare = dipolePairDensityFromIcSquareToRoiSquare(currentSessionIdx, :);
    
    %Get normalization factor and normalize
    normFactor = sum(dipolePairDensitiesInIcSquare);
    dipolePairDensitiesInIcSquareNormalized = dipolePairDensitiesInIcSquare ./ normFactor;
    
    %Replace NaNs with zeros caused by division by zero
    dipolePairDensitiesInIcSquareNormalized(isnan(dipolePairDensitiesInIcSquareNormalized)) = 0;

    effectiveConnectivityTimeFreq = single(dipolePairAndMeasureObj.linearizedMeasure(currentSessionIdx,:)' * dipolePairDensitiesInIcSquareNormalized);
    
    %Reshape and re-order data
    effectiveConnectivityTimeFreq = reshape(effectiveConnectivityTimeFreq,  timeFreqSize(1), timeFreqSize(2), size(effectiveConnectivityTimeFreq,2));
    effectiveConnectivityTimeFreq = permute(reshape(effectiveConnectivityTimeFreq, timeFreqSize(1), timeFreqSize(2), max(fromRegionNumber), max(fromRegionNumber)), [3 4 1 2]); % Use most consistent dimensions for SIFT and most intuitive for users
  
    countMatrix(:,:,sessionIdx) = logical(squeeze(sum(sum(effectiveConnectivityTimeFreq,3),4)));
    disp(horzcat('Checked file ', num2str(sessionIdx)))
end
sumEdgeMatrix        = sum(countMatrix,3);
numSubjThreshold     = size(countMatrix,3)*userInputPercentage;
stillGoodMask        = sumEdgeMatrix>=numSubjThreshold;


%Now do actual connectivity projection, subj thresholding, and saving
for sessionIdx = unique(dipolePairAndMeasureObj.sessionId)'
    tic
    %Will use matrix multiplication instead of old way
    currentSessionIdx = find(dipolePairAndMeasureObj.sessionId==sessionIdx); 
    dipolePairDensitiesInIcSquare = dipolePairDensityFromIcSquareToRoiSquare(currentSessionIdx, :);
    
    %Get normalization factor and normalize
    normFactor = sum(dipolePairDensitiesInIcSquare);
    dipolePairDensitiesInIcSquareNormalized = dipolePairDensitiesInIcSquare ./ normFactor;
    
    %Replace NaNs with zeros caused by division by zero
    dipolePairDensitiesInIcSquareNormalized(isnan(dipolePairDensitiesInIcSquareNormalized)) = 0;

    effectiveConnectivityTimeFreq = single(dipolePairAndMeasureObj.linearizedMeasure(currentSessionIdx,:)' * dipolePairDensitiesInIcSquareNormalized);
    
    %Reshape and re-order data
    effectiveConnectivityTimeFreq = reshape(effectiveConnectivityTimeFreq,  timeFreqSize(1), timeFreqSize(2), size(effectiveConnectivityTimeFreq,2));
    effectiveConnectivityTimeFreq = permute(reshape(effectiveConnectivityTimeFreq, timeFreqSize(1), timeFreqSize(2), max(fromRegionNumber), max(fromRegionNumber)), [3 4 1 2]); % Use most consistent dimensions for SIFT and most intuitive for users

    %Cut out edges that don't need subject threshold criteria
    %effectiveConnectivityTimeFreq = bsxfun(@times, effectiveConnectivityTimeFreq, groupEdgesMeetingThreshold);
    effectiveConnectivityTimeFreq = bsxfun(@times, effectiveConnectivityTimeFreq, stillGoodMask);
    
    %Store only non-zero values to save space
    nonzeroConn = effectiveConnectivityTimeFreq(effectiveConnectivityTimeFreq~=0);
    
    %Since if a connection exists, data will be available for any
    %time/freq, then only need to store idx of first time/freq 
    subjectEdgeIdx = find(sum(sum(effectiveConnectivityTimeFreq,3),4)~=0); %nonzeroIdx for 
    
    %finallySelectedEdgeIdx   = find(edgesMeetingThreshold);
    dims = size(effectiveConnectivityTimeFreq);
    if length(dims) == 3
        dims(4) = 1;
    end

    %This fcn will get the full list of idx of the stored data
        %nonzeroIdxFull = reshape(bsxfun(@plus, double(finallySelectedEdgeIdx), dims(1)*dims(2)*[0:dims(3)*dims(4)-1]),[],1);
    %Then will re-create the original large matrix (n x n x freq x time)
        %connMatFull = zeros(dims);
        %connMatFull(nonzeroIdxFull) = nonzeroConn;
    
    %Save individual subject files
    connectivityType = 'rPDC';
    originalFile = fileNameList{sessionIdx};
    parcelInfo = networkTxt;
    finallySelectedEdgeIdx = find(stillGoodMask);
    save([workingFolder subIDList{sessionIdx} '.conn'],...
        'nonzeroConn', 'subjectEdgeIdx', 'finallySelectedEdgeIdx',...
        'dimensionLabels', 'dims', 'frequencies', 'latencies',...
        'connectivityType', 'fileNameList',...
        'parcelInfo', 'originalFile', 'userInputPercentage', '-v6')
    
    timeLapse = toc;
    disp(sprintf('%2.0d/%2.0d subjects done (%0.1d sec lapsed)', sessionIdx, length(unique(dipolePairAndMeasureObj.sessionId)), round(timeLapse)));

end
1;


% 12/09/2019 Makoto. When dual dipoles are detected, compare their moments and pick up the larger one.
function output = selectDipWithLargerMoment(posxyz, momxyz)

if norm(posxyz)> 200
    error('Abnormal dipole location detected.')
end

if size(posxyz, 1) == 2
    mom1 = norm(momxyz(1,:));
    mom2 = norm(momxyz(2,:));
    if mom1 > mom2
        output = posxyz(1,:);
    else
        output = posxyz(2,:);
    end
else
    output = posxyz;
end

end