% This script runs wide-band connectivity analysis of 3 methods: xPAC, PLV
% and Theta Power.
% wide-band connections of theta and gamma are slow-slow(2-5 to 30-70),
% slow-fast (2-5 to 70-100),fast-slow (5-9 to 30-70), and fast-fast(30-70 to 30-70)
% xPAC by Canolty with 250 surrogate data;
% PLV by Laux.
% theta power buy normalized theta power,(faraction of theta power over full-band spectral power)
%
% necessary functions are fund in path folders whereas iEEG data is stored
% in EEG_DATA, pre-processed by Sarah.
%
%
% Functions are desgined to running the data locally
% Fuctions used : ROIabbr, xPACMI_ele, whohasROI, and ttest_xPAC_ele.
%
% David Wang
% latest update: 11/2/2019

clear
clc;close all;



% all subjects has talStruct

subjects_total =  {'UT021','UT027','UT049','UT050','UT053','UT056','UT063','UT064','UT067','UT069',...
    'UT074','UT083','UT084','UT090','UT095','UT113','UT114','UT117','UT122','CC013','CC014','CC015',...
    'CC016','CC017','UT008','UT009','UT011','UT013','UT018','UT020','UT021','UT025','UT034',...
    'UT037','UT039','UT077','UT079','UT081','UT103','UT111'};

sigacc = 1;
insigacc=1;
phase_su = cell(1,1);
phase_insig = cell(1,1);
sig_f = cell(1,1);
f_insig = cell(1,1);
meanMI = cell(1,1);
ROIs_t = {'AH','PH','EC','paraHG','LPF','LMT','PC','BTL','SM'};
hemit = {'_L','_R'}; % hemisphere chosen from _L and _R



for Rind = 1:length(ROIs_t)
    
    
    ThisROI  = ROIs_t(Rind);
    
    tic
    for a = 1:length(hemit)
        eleacc =1;
        hemi = hemit(a);
        ROIs =ROIabbr(string(strcat(ThisROI,hemi))); % left hemisphere
        
        for subind = 1:length(subjects_total)
            subjects = subjects_total(subind);
            
            load info.mat
            load ROIabbre.mat
            Nochannel = whohasROI(subjects,ROIs);
            
            Fs = 250;
            
            load(string(subjects));
            
            index_recalled = find([events.recalled]==1);
            index_nonrecalled = find([events.recalled]==0);
            recalled_eeg = eeg(:,index_recalled,:);
            nonrecalled_eeg = eeg(:,index_nonrecalled,:);
            
            indx =  str2num(Nochannel{1});
            ind_total = find(ismember(channels, indx));
            
            if isempty(Nochannel{1})
                continue
            end
            
            
            
            for base_ind = 1:length(ind_total)
                ThisRegion{1} = reshape(recalled_eeg(ind_total(base_ind),:,:),[size(recalled_eeg,2),size(recalled_eeg,3)]);
                ThisRegion{2} = reshape(nonrecalled_eeg(ind_total(base_ind),:,:),[size(nonrecalled_eeg,2),size(nonrecalled_eeg,3)]);
                
                
                df = [3 4];
                df2 = 40;
                f1 = [2 5];
                f2 = [30 70];
                
                for j = 1:length(f2)
                    band = [f2(j),f2(j)+df2];
                    phase_theta = cell(1,2);
                    P_norm = cell(1,2)';
                    for cond = 1:2
                        band = [f2(j),f2(j)+df2];                       
                        ThisCond = ThisRegion{cond};
                        power_all = pwer(ThisCond);
                        power_gamma = pwer(PWprep(ThisCond,band,Fs));
                        
                        P_norm{cond} = power_gamma./power_all;
                    end
                    
                    
                    SME_su{eleacc}{1,j} = P_norm{1};
                    SME_unsu{eleacc}{1,j} = P_norm{2};
                end
                
                
                eleacc = eleacc+1;
                
                
            end
            
        end
        
        path2save = string(strcat(cd,'/wb/SME/'));
        if ~exist(path2save,'dir')
            mkdir(path2save)
        end
        
        whatRcomputing = string(strcat(base,'_',ROI,hemi));
        sprintf(whatRcomputing)
        
        path2save1 = string(strcat(cd,'/wd/SME/'));
        file2save1 = strcat(path2save1,whatRcomputing,'.mat');
        if ~exist(path2save1,'dir')
            mkdir(path2save1)
        end
        
        SME_total = [SME_su;SME_unsu];
        save(string(file2save1),'SME_total');
        clear SME_su SME_unsu
        toc
        
    end
    
end


sprintf('finished')