% wide-band results processing, circular connection plots are generated in
% plotconn.m script


% data matrix cell form: [slow-slow,slow-fast;fast-slow,fast-fast]

clc
clear all
close all


addpath('/Volumes/Study/PAC/statistics');
addpath('/Volumes/Study/PAC/wide-band/wb/PLV')
addpath('/Volumes/Study/PAC/wide-band/wb/xPAC')
addpath('/Volumes/Study/PAC/plotting_tools/circro-master')



baset= {'AH','PH'};
hemit = {'_L','_R'};
ROIs = {'EC','paraHG','LPF','LMT','PC','BTL','SM'};


labelL = strcat([baset,ROIs],hemit(1));
labelR = strcat([baset,ROIs],hemit(2));
Labels = [labelL',labelR'];


cal = 'xPAC';
if strcmp(cal,'xPAC')
    band = {'slow','slowfast','fastslow','fast'};
else
    band = {'slow','fast'};
end



ratio_matrix= cell(1,2);
ratio_matrix{1} = zeros(9,9);
ratio_matrix{2} = zeros(9,9);
myratio = cell(1,length(band));


NoElePair = zeros(length(ROIs),4);
acc= 1;


for baseind = 1:length(baset)
    
    base = baset(baseind);
    
    for hind = 1:length(hemit)
        hemi = hemit(hind);
        
        
        for Rind = 1:length(ROIs)
            
            
            
            
            ROI = ROIs(Rind);
            
            whatRretriving = string(strcat('/Volumes/Study/PAC/wide-band/wb/',cal,'/',base,'_',ROI,hemi,'.mat'));
            
            load(whatRretriving)
            
            
            NoElePair(Rind,acc) = length(MI_total);                    
        end
        acc= acc+1;
    end        
end

NoElePair

minNo = min(NoElePair(:))
maxNo = max(NoElePair(:))
