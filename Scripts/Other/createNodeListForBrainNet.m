%Create .node file for BrainNet for gordon mapping
[~,~,networkTxt] = xlsread('D:\NetworkConnectivity\Gordon2016_Map\Parcels\Parcels.xlsx');

networkNames = unique(networkTxt(:,6));

%Exclude any networks (e.g., excludeNetworks = [5])
excludeNetworks               = [];
networkNames(excludeNetworks) = [];

newParcelOrder = [];
for i = 1:size(networkNames,1)
    parcelsIdxNetwork = find(strcmp(networkTxt(:,6),networkNames(i))==1);
    newParcelOrder = [newParcelOrder ; parcelsIdxNetwork];
end

newNetwork = networkTxt(newParcelOrder,:);

nodeArray = {};

for i=1:size(newNetwork,1)
    nodeArray{i,1} = horzcat(newNetwork{i,5},' 1 1 ', newNetwork{i,6}, newNetwork{i,3});
end

%Save nodeArray to .node file