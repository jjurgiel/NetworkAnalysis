[group1Files, group1Folder] = uigetfile('*.conn', 'MultiSelect', 'on');
if ~any(group1Folder); disp('Cancelled.'); return; end

%Select group 2 data
[group2Files, group2Folder] = uigetfile('*.conn', 'MultiSelect', 'on');
if ~any(group1Folder); disp('Cancelled.'); return; end

group1Array = [];
group2Array = [];

for n = 1:size(group1Files, 2)
    load(horzcat(group1Folder, group1Files{n}), '-mat');
    group1Array(n) = length(subjectEdgeIdx);
end

for n = 1:size(group2Files, 2)
    load(horzcat(group2Folder, group2Files{n}), '-mat');
    group2Array(n) = length(subjectEdgeIdx);
end