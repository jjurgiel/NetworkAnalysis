% This function takes in a AAL region and finds the Gordon parcel which 
% closest corresponds to it.
%
% This is necessary as the Gordon parcellation does not have nice
% anatomical labels.
%
% Gordon Parcellations:
% Gordon, Evan M., et al. "Generation and evaluation of a cortical area 
% parcellation from resting-state correlations." Cerebral cortex 26.1 (2016): 288-303.
%
% ROIs available (alphabetically listed):
% 'Angular.L';              'Angular.R';            'Calcarine.L';
% 'Calcarine.R';            'Cingulum.Ant.L';       'Cingulum.Ant.R';
% 'Cingulum.Mid.L';         'Cingulum.Mid.R';       'Cingulum.Post.L';
% 'Cingulum.Post.R';        'Cuneus.L';             'Cuneus.R';
% 'Frontal.Inf.Oper.L';     'Frontal.Inf.Oper.R';   'Frontal.Inf.Orb.L';
% 'Frontal.Inf.Orb.R';      'Frontal.Inf.Tri.L';    'Frontal.Inf.Tri.R';
% 'Frontal.Med.Orb.L';      'Frontal.Med.Orb.R';    'Frontal.Mid.L';
% 'Frontal.Mid.Orb.L';      'Frontal.Mid.Orb.R';    'Frontal.Mid.R';
% 'Frontal.Sup.L';          'Frontal.Sup.Medial.L'; 'Frontal.Sup.Medial.L';
% 'Frontal.Sup.Orb.L';      'Frontal.Sup.Orb.R';    'Frontal.Sup.R';
% 'Fusiform.L';             'Fusiform.R';           'Insula.L';
% 'Insula.R';               'Lingual.L';            'Lingual.R';
% 'LowerBasal.L';           'LowerBasal.R';         'Occipital.Inf.L';
% 'Occipital.Inf.R';        'Occipital.Mid.L';      'Occipital.Mid.R';
% 'Occipital.Sup.L';        'Occipital.Sup.R';      'Paracentral.Lobule.L';
% 'Paracentral.Lobule.R';   'Parietal.Inf.L';       'Parietal.Inf.R';
% 'Parietal.Sup.L';         'Parietal.Sup.R';       'Postcentral.L';
% 'Postcentral.R';          'Precentral.L';         'Precentral.R';
% 'Precuneus.L';            'Precuneus.R';          'Rectus.L';
% 'Rectus.R';               'Rolandic.Oper.L';      'Rolandic.Oper.R';
% 'Supp.Motor.Area.L';      'Supp.Motor.Area.R';    'SupraMarginal.L';
% 'SupraMarginal.R';        'Temporal.Inf.L';       'Temporal.Inf.R';
% 'Temporal.Mid.L';         'Temporal.Mid.R';       'Temporal.Pole.Mid.L';
% 'Temporal.Pole.Mid.R';    'Temporal.Pole.Sup.L';  'Temporal.Pole.Sup.R';
% 'Temporal.Sup.L';         'Temporal.Sup.R';       'UpperBasal.L';
% 'UpperBasal.R'
%
%
% Written by Joseph Jurgiel @ UCLA Semel Institute (jjurgiel@ucla.edu)
% July 2022

%Select ROI
selectedROI = 'Frontal.Sup.R';

%%
%Load in .node files containing coordinates
tableAAL = table2cell(readtable('D:\NetworkConnectivity\AAL3v1_for_SPM12\AAL90jj.node', 'filetype', 'text'));
tableGordon = table2cell(readtable('D:\NetworkConnectivity\Gordon2016_Map\Gordon286Nodes.node', 'filetype', 'text'));

%Find idx of ROI
[idxROI, ~] = find(strcmp(selectedROI, tableAAL));

%Find coordinates of ROI, and coordinates of parcels
coordAAL = cell2mat(tableAAL(idxROI,1:3));
coordGordon = cell2mat(tableGordon(:,1:3));

%Find parcel closest to ROI
dists = vecnorm(coordGordon - coordAAL,2,2);
parcelIdx = find(dists == min(dists));

disp(['Closest parcel for ', selectedROI, ' is Parcel #', num2str(parcelIdx), '(', num2str(min(dists)), '): ', tableGordon{parcelIdx,6}]);