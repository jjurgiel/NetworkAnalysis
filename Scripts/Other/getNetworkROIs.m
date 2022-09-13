% This compares overlap between nii parcelated file and groupSIFT ROI grid
% Jan 2022, for Giorgia & Sandy theorizing network connectivity
%{
if ~exist('dOut','var')
    ALLEEG = [];
    EEG = pop_loadset('filename', '5001 pruned with ICAICremovedStim-Check-CongruCorrect.set',...
                      'filepath','D:\forJoe\Project1\TestData\Original\CongruCorrect', 'loadmode', 'info');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );  


    dOut = getGridLocations(ALLEEG,50);
end


roiGridLocs = []; 
for i = 1:size(dOut.roiLabels,1)
    roiOut =  pr.regionOfInterestFromAnatomy(pr.headGrid(2), dOut.roiLabels{i});

    roiIdx = find(roiOut.membershipCube==1);
    
    gridLocsX = roiOut.headGrid.xCube(roiIdx);
    gridLocsY = roiOut.headGrid.yCube(roiIdx);
    gridLocsZ = roiOut.headGrid.zCube(roiIdx);
    
    roiGridLocs = [roiGridLocs;[repmat(i,size(roiIdx,1),1) gridLocsX gridLocsY+1 gridLocsZ+1]];
    disp(horzcat('Finished ROI',num2str(i)))
end
%}

%niiMatrix = niftiread('C:\Users\jajur\Downloads\Parcels-19cwpgu\Parcels\Parcels_MNI_222.nii');
niiInfo = niftiinfo( 'D:\NetworkConnectivity\1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels\MNI152\Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz');
niiMatrix = niftiread('D:\NetworkConnectivity\1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels\MNI152\Yeo2011_7Networks_N1000.split_components.FSL_MNI152_2mm.nii.gz');
xRange = -90:2:90; yRange = -126:2:90; zRange = -72:2:108;

%Get parcel of each grid point
for n = 1:size(roiGridLocs,1)
    xIdx = find(xRange == roiGridLocs(n,2));
    yIdx = find(yRange == roiGridLocs(n,3));
    zIdx = find(zRange == roiGridLocs(n,4));
    
    roiGridLocs(n,5) = niiMatrix(xIdx,yIdx,zIdx);
    
    if mod(n,2000) == 0
       disp(horzcat('Finished grid point #',num2str(n)))
    end
end

%Check which parcels overlap each ROI
parcelList = {};
for r = 1:76
    gridIdx = find(roiGridLocs(:,1)==r);
    subsetGrid = roiGridLocs(gridIdx,:);
    
    uniqueParcels = unique(subsetGrid(:,5));
    for n = 1:length(uniqueParcels)
        numParcelPts = length(find(subsetGrid(:,5)==uniqueParcels(n)));
        if uniqueParcels(n) ~= 0
            parcelList{r}{n,1} = dOut.roiLabels{uniqueParcels(n)};
        else
            parcelList{r}{n,1} = 'None';
        end
        parcelList{r}{n,2} = numParcelPts;
    end
end

ROIList = {};
%Check which ROIs overlap each parcel
for p = 1:length(unique(roiGridLocs(:,5)))-1
    gridIdx = find(roiGridLocs(:,5)==p);
    subsetGrid = roiGridLocs(gridIdx,:);
    
    %Get total # of grid points in parcel
    totalParcelGridCount(p) = sum(sum(sum(niiMatrix==p)));
    
    %Make list of ROI %s in each parcel
    uniqueROIs = unique(subsetGrid(:,1));
    for n = 1:length(uniqueROIs)
        numROIPts = length(find(subsetGrid(:,1)==uniqueROIs(n)));
        ROIList{n,p} = horzcat(dOut.roiLabels{uniqueROIs(n)}, ' ', num2str(round(100*numROIPts/totalParcelGridCount(p))),'%');
      %  ROIList{p}{n,2} = numROIPts;
    end
    
    
end

%Check number of grid points with no parcel label
noneArray = [0, 11, 18, 19, 73, 115, 118:125, 128, 129, 133:135, 142, 144, 172, 178, 179, 280:289, 291, 292, 296, 297, 300:306,312,314,];
noneSum = 0;
for n = 1:size(roiGridLocs,1)
    if ismember(roiGridLocs(n,5),noneArray)
       noneSum = noneSum+1; 
    end
end