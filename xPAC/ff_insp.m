%frequency (Hippocampus) - frquency(ROI) inspection

% This is the main script of Cross-region Phase-Amplitude Coupling(xPAC)
% of subjects and ROIs at electrodes level.
% differ from main_xPAC.m file, this script adapts core functions of
% xPACMI_ele and ttest_xPAX_ele at electrodes level as well.
%
%
% Functions are desgined to running the data locally
% Fuctions used : ROIabbr, xPACMI_ele, whohasROI, and ttest_xPAC_ele.

% David Wang 03/20/2019
clear all
clc;close all; 

addpath('/Volumes/Study/EEG_DATA/');
addpath('/Volumes/Study/PAC/');
addpath('/Volumes/Study/PAC/PLV')

% all subjects has talStruct
subjects =  {'UT027'};



load info.mat
load ROIabbre.mat

hemi = '_L'; % hemisphere chosen from _L and _R
base = 'AH'; % base region that theta phase comes from, chosen from AH or PH
band = 'slowfast';
other_ROIs = {'paraHG'};


ROIs  = ['AH', other_ROIs];



ROIs =ROIabbr(strcat(ROIs,hemi)); % left hemisphere

[subjhas, Nochannel] = whohasROI(subjects,ROIs);

Fs = 250;

load UT027;

index_recalled = find([events.recalled]==1);
index_nonrecalled = find([events.recalled]==0);
recalled_eeg = eeg(:,index_recalled,:);
nonrecalled_eeg = eeg(:,index_nonrecalled,:);



ah(:,1,:) = (recalled_eeg(:,11,:));
ah(:,2,:) = (nonrecalled_eeg(:,11,:));

para(:,1,:) = (recalled_eeg(:,10,:));
para(:,2,:) = (nonrecalled_eeg(:,10,:));

tic
f1 = 2:1:7;
f2 = 40:1:99;

for i = 1:length(f1)
    for j = 1:length(f2)
        
        band = [f1(i)-1, f1(i)+1; f2(j)-1,f2(j)+1];
        amp_gamma = zeros(126,450);
        phase_theta = zeros(126,450);
        M_raw=cell(1,2);
        
        for cond = 1:2
            AH = reshape(ah(:,cond,:),[126,450]);
            PARA = reshape(para(:,cond,:),[126,450]);
            
            for k = 1:126
                [amp_gamma(i,:),~]=PACprep(AH(i,:),band,Fs);
                [~,phase_theta(i,:)]=PACprep(PARA(i,:),band,Fs);
            end
            
            M_raw{cond} = PLV(amp_gamma,phase_theta);
            pac{cond}(i,j)=mean(M_raw{cond});              
        end
        [~,pvalue(i,j)] = ttest2(M_raw{1},M_raw{2},'tail','both');
    end    
end
toc


figure
subplot(1,3,1)
contour(f1,f2,pac{1}')
set(gca,'YDir','normal')
colorbar
subplot(1,3,2)
contour(f1,f2,pac{2}')
set(gca,'YDir','normal')
colorbar
subplot(1,3,3)
% 
% 
% pcolor(f2,f1,pvalue)
% set(gca,'YDir','normal')
%hold on
contour(f1,f2,pvalue')
set(gca,'YDir','normal')
ylabel('paraHG');xlabel('AH')
colorbar
title('pvalue of AH-paraHG between successful and unsuccessful encoding')

