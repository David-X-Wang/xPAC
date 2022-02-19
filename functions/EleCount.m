
clc;close all; clear

addpath('/Volumes/David 5TB/EEG_DATA/');
addpath('/Volumes/Study/PAC/');
addpath('/Volumes/Study/PAC/plotting_tools/cbrewer')
% all subjects has talStruct 
subjects =  {'UT021','UT027','UT049','UT050','UT053','UT056','UT063','UT064','UT067','UT069',...
    'UT074','UT083','UT084','UT090','UT095','UT113','UT114','UT117','UT122','CC013','CC014','CC015',...
    'CC016','CC017','UT008','UT009','UT011','UT013','UT018','UT020','UT021','UT025','UT034',...
    'UT037','UT039','UT077','UT079','UT081','UT103','UT111'};


load info.mat
load ROIabbre.mat

hemi_t = {'_L','_R'};

eleNoforRegion = zeros(9,2);

for k = 1:length(hemi_t)
  hemi = hemi_t(k);   
%hemi = '_L'; % hemisphere chosen from _L and _R
base = 'AH'; % base region that theta phase comes from, chosen from AH or PH
ROIs = {'AH','PH','EC','paraHG','LPF','LMT','PC','BTL','SM'};
% 
% if strcmp(base,'AH')
%     ROIs  = ['AH','PH', other_ROIs];
% else
%     ROIs  = ['PH','AH', other_ROIs];
% end


ROIs =ROIabbr(strcat(ROIs,hemi)); % left hemisphere

[Nochannel] = whohasROI(subjects,ROIs);
eleNo = zeros(size(Nochannel));
for i = 1:size(Nochannel,1)
    for j = 1:size(Nochannel,2)
      eleNo(i,j) = length(cell2mat(Nochannel(i,j)));
    end
end

eleNoforRegion(:,k) = sum(eleNo')';
end
n = (max(eleNoforRegion(:)));
cmap = (cbrewer('seq','YlOrRd',n)); % blues at bottom

figure('position',[0 0 500 500])
image(eleNoforRegion)
colormap(cmap)
 caxis([min(eleNoforRegion(:)) max(eleNoforRegion(:))])
h = colorbar('position',[0.92 0.41 0.03 0.2]);
% print('-painters','-depsc','/Volumes/Study/PAC/Paper/figures/Fig2mat/colorbar.eps')
caxis([min(eleNoforRegion(:)) max(eleNoforRegion(:))])
    
colorsL= round(255*cmap(eleNoforRegion(:,1),:));
colorsR = round(255*cmap(eleNoforRegion(:,2),:));

