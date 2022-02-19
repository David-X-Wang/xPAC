% This script runs wide-band connectivity computation of xPAC and PLV.
% For xPAC (by Canolty with 250 surrogate data):
% theta and gamma are split into slow-slow(2-5 Hzto 30-70 Hz),slow-fast (2-5 to 70-100Hz)
% fast-slow (5-9 to 30-70Hz), and fast-fast(30-70 to 30-70Hz)
%
% For PLV (by Laux.):
% theta is split into slow (2-5 to 2-5 Hz) and fast (5-9Hz to 5-9Hz) bands.
%
% Precursor functions are found in path folders
% iEEG data (channels/regions*trials/tasks*samples/time)is stored in
% EEG_DATA folder, pre-processed by Sarah Seger.
%
%
% Functions are desgined to running the data locally
% Fuctions used : ROIabbr, xPACMI_ele, whohasROI, and ttest_xPAC_ele.
%
% David Wang
% latest update: 11/5/2019

clear
clc;close all;


% all subjects has talStruct

subjects_total =  {'UT021','UT027','UT049','UT050','UT053',...
    'UT056','UT063','UT064','UT067','UT069',...
    'UT074','UT083','UT084','UT090','UT095',...
    'UT113','UT114','UT117','UT122'};


sigacc = 1; % accumulator for significant  pairs
insigacc=1; % accumulator for nonsiginificant pairs
phase_su = cell(1,1);
phase_insig = cell(1,1);
sig_f = cell(1,1);
f_insig = cell(1,1);
meanMI = cell(1,1);
ROIs_total = {'paraHG'};%,'paraHG','LPF','LMT','PC','BTL','SM'};
hemi_total = {'_L','_R'};
base_total = {'AH','PH'};


for Rind = 1:length(ROIs_total)
    ROI = ROIs_total(Rind);
    tic
    for a = 1:length(hemi_total)
        
        for b= 1:length(base_total)
            eleacc =1; % accumulator of electrode pairs for the connection
            for subind = 1:length(subjects_total)
                subjects = subjects_total(subind);
                load info.mat
                load ROIabbre.mat
                hemi = hemi_total(a);
                base = base_total(b);
                % form hippocampus-ROI connection
                ROIs  = [base, ROI];
                % identify full name of the region in the hemisphere. (e.g. AH-L to Anterior Hippocampus Left)
                ROIs =ROIabbr(strcat(ROIs,hemi));
                % extract channel numbers of the regions for the subjects.
                Nochannel = whohasROI(subjects,ROIs);
                
                Fs = 250;
                
                load(string(subjects));
                
                % Seperate recalled (SuE) and nonrecalled (unsuE) form the
                % encoding EEGs.
                index_recalled = find([events.recalled]==1);
                index_nonrecalled = find([events.recalled]==0);
                recalled_eeg = eeg(:,index_recalled,:);
                nonrecalled_eeg = eeg(:,index_nonrecalled,:);
                
                base_indx =  str2num(Nochannel{1});
                ROI_indx = str2num(Nochannel{2});
                
                channels = [talStruct.channel]';
                
                base_ind_total = find(ismember(channels, base_indx));
                ROI_ind_total=find(ismember(channels,ROI_indx));
                % skip to next subject if it does not have the base (AH/PH)
                if isempty(Nochannel{1})
                    continue
                end
                
                % skip to next subject if it does not have the ROI
                if isempty(Nochannel{2})
                    continue
                end
                
                % run the xPAC of all possible electrode pairs (combinations of two regions)
                
                for base_ind = 1:length(base_ind_total)
                    for ROI_ind = 1:length(ROI_ind_total)
                        
                        % reshape the base region trials of the electrode
                        % to the form : trials*samples, {1} for recaleld,
                        % {2} for nonrecalled
                        ThisBase{1} = reshape(recalled_eeg(base_ind_total(base_ind),:,:),[size(recalled_eeg,2),size(recalled_eeg,3)]);
                        ThisBase{2} = reshape(nonrecalled_eeg(base_ind_total(base_ind),:,:),[size(nonrecalled_eeg,2),size(nonrecalled_eeg,3)]);
                        
                        % reshape the base region trials of the electrode
                        % to the form : trials*samples
                        ThisROI{1} = reshape(recalled_eeg(ROI_ind_total(ROI_ind),:,:),[size(recalled_eeg,2),size(recalled_eeg,3)]);
                        ThisROI{2} = reshape(nonrecalled_eeg(ROI_ind_total(ROI_ind),:,:),[size(nonrecalled_eeg,2),size(nonrecalled_eeg,3)]);
                        
                        
                        df = 5;
                        df2 = 40;
                        f1 = 2:df:11;
                        f2 = 30:df2:100;
                        pac_norm = cell(1,6);
                        
                        for i = 1:length(f1)
                            
                            band = [f1(i), f1(i)+df];
                            amp_gamma = cell(1,2);
                            phase_theta = cell(1,2);
                            
                            phase_theta_ROIs = cell(1,2);
                            M_norm = cell(1,2);
                            pfPha = cell(1,2);
                            PLVs = cell(1,2);
                            
                            
                            for cond = 1:2
                                ThisConBase = ThisBase{cond};
                                ThisConROI = ThisROI{cond};
                                phase_theta=PLVprep(ThisConBase,band,Fs); % theta phases for Base               
                                phase_theta_ROIs=PLVprep(ThisConROI,band,Fs);
       
                                % PLV
                                PLVs{cond} = PLV(phase_theta,phase_theta_ROIs);
                            end
                            
                            PLV_su{eleacc}{i} = PLVs{1};
                            PLV_unsu{eleacc}{i} = PLVs{2};
                            
                        end
                        eleacc = eleacc+1;
                    end
                end
                
            end
            whatRcomputing = string(strcat(base,'_',ROI,hemi));
            sprintf(whatRcomputing)
            
            path2save1 = string(strcat(cd,'/wd/PLV/'));
            file2save1 = strcat(path2save1,whatRcomputing,'.mat');
            if ~exist(path2save1,'dir')
                mkdir(path2save1)
            end
            
            PLV_total = [PLV_su;PLV_unsu];
            
            
            save(string(file2save1),'PLV_total');
            
            
            clear xPAC_p PLV_su PLV_unsu
            toc
        end
    end
end


sprintf('finished')