%frequency (Hippocampus) - frquency(ROI) inspection
% xPAC by Canolty, exclusively
% This is the main script of frequency to frequency inspection of (xPAC)
% of subjects and ROIs at electrodes level.
% differ from main_xPAC.m file, this script adapts core functions of
% xPACMI_ele and ttest_xPAX_ele at electrodes level as well.
%
%
% Functions are desgined to running the data locally
% Fuctions used : ROIabbr, xPACMI_ele, whohasROI, and ttest_xPAC_ele.

% David Wang 10/18/2019
clear
clc;close all;
tic
if ismac
    addpath('/Volumes/Study/EEG_DATA/');
    addpath('/Volumes/Study/PAC/xPAC');
    addpath('/Volumes/Study/PAC/PLV');
    addpath('/Volumes/Study/PAC/Power');
    addpath('/Volumes/Study/PAC/');
elseif isunix
    addpath('/users/xiaoliangw/EEG_DATA');
    addpath('/users/xiaoliangw/PAC/statistics');
    addpath('/users/xiaoliangw/PAC/xPAC');
    addpath('/users/xiaoliangw/PAC/PLV');
    addpath('/users/xiaoliangw/PAC/Power');
    addpath('/users/xiaoliangw/PAC/');
end

% all subjects has talStruct

subjects_total =  {'UT021','UT027','UT049','UT050','UT053',...
    'UT056','UT063','UT064','UT067','UT069',...
    'UT074','UT083','UT084','UT090','UT095',...
    'UT113','UT114','UT117','UT122'};
sigacc = 1;
insigacc=1;
phase_su = cell(1,1);
phase_insig = cell(1,1);
sig_f = cell(1,1);
f_insig = cell(1,1);
meanMI = cell(1,1);

for subind = 1:length(subjects_total)
    subjects = subjects_total(subind);
    
    cal = 'xPAC';
    
    load info.mat
    load ROIabbre.mat
    
    hemit = {'_L','_R'}; % hemisphere chosen from _L and _R
    baset = {'AH','PH'}; % base region that theta phase comes from, chosen from AH or PH
    
    for a = 1:length(hemit)
        for b= 1:length(baset)
            hemi = hemit(a);
            base = baset(b);
            
            other_ROIs = {'paraHG'};
            ROIs  = [base, other_ROIs];
            ROIs =ROIabbr(strcat(ROIs,hemi)); % left hemisphere
            
            Nochannel = whohasROI(subjects,ROIs);
            
            Fs = 250;
            
            load(string(subjects));
            
            index_recalled = find([events.recalled]==1);
            index_nonrecalled = find([events.recalled]==0);
            recalled_eeg = eeg(:,index_recalled,:);
            nonrecalled_eeg = eeg(:,index_nonrecalled,:);
            
            base_indx =  str2num(Nochannel{1});
            ROI_indx = str2num(Nochannel{2});
            
            base_ind_total = find(ismember(channels, base_indx));
            ROI_ind_total=find(ismember(channels,ROI_indx));
            
            paircount = 1;
            if isempty(Nochannel{1})
                continue
            end
            if isempty(Nochannel{2})
                continue
            end
            
            for base_ind = 1:length(base_ind_total)
                for ROI_ind = 1:length(ROI_ind_total)
                    
                    
                    
                    ah{1} = reshape(recalled_eeg(base_ind_total(base_ind),:,:),[size(recalled_eeg,2),size(recalled_eeg,3)]);
                    ah{2} = reshape(nonrecalled_eeg(base_ind_total(base_ind),:,:),[size(nonrecalled_eeg,2),size(nonrecalled_eeg,3)]);
                    
                    para{1} = reshape(recalled_eeg(ROI_ind_total(ROI_ind),:,:),[size(recalled_eeg,2),size(recalled_eeg,3)]);
                    para{2} = reshape(nonrecalled_eeg(ROI_ind_total(ROI_ind),:,:),[size(nonrecalled_eeg,2),size(nonrecalled_eeg,3)]);
                    
                    tic
                    df = 1 ;
                    df2 = 2;
                    f1 = 2:df:11;
                    f2 = 40:df2:100;
                     pac_norm = cell(1,5);
                    for i = 1:length(f1)
                        for j = 1:length(f2)
                            
                            
                            band = [f1(i), f1(i)+df; f2(j),f2(j)+df2];
                            
                            %M_raw=cell(1,2);
                            
                            amp_gamma = cell(1,2);
                            
                            phase_theta = cell(1,2);
                            phase_theta_ROIs = cell(1,2);
                            M_norm = cell(1,2);
                            pfPha = cell(1,2);
                            PLVs = cell(1,2);
                            P_norm = cell(1,2);
                           
                                
                            for cond = 1:2
                                AH = ah{cond};
                                PARA = para{cond};
                                

                                %MInorm = zeros(1,size(ah{cond},1));
                                
                                % PAC 
                                [amp_gamma{cond},phase_theta_ROIs{cond}]=PACprep(PARA,band,Fs);
                                [~,phase_theta{cond}]=PACprep(AH,band,Fs); 
                                [M_norm{cond}, pfPha{cond}] = surr_max(amp_gamma{cond},phase_theta{cond});
                           
                                % PLV
                                PLVs{cond} = PLV(phase_theta{cond},phase_theta_ROIs{cond});
                                
                                % Power
                                power_all = pwer(AH);
                                power_theta = pwer(PWprep(AH,band,Fs));
                                P_norm{cond} = power_theta./power_all;
                                
                                
                                pac_norm{cond}(i,j)=mean(M_norm{cond});


                            end
                            
                            % T-test of PAC, PLV and Power connectivity
                            % between SuE and UnsuE.
                            [~,pvalue1] = ttest2(M_norm{1}(isfinite(M_norm{1})),M_norm{2}(isfinite(M_norm{2})),'tail','both');
                            [~,pvalue2] = ttest2(PLVs{1}(isfinite(PLVs{1})),PLVs{2}(isfinite(PLVs{2})),'tail','both');                            
                            [~,pvalue3] = ttest2(P_norm{1}(isfinite(P_norm{1})),P_norm{2}(isfinite(P_norm{2})),'tail','both');
                            % find those significant electrode pairs
                            pac_norm{3}(i,j) = pvalue1;
                            pac_norm{4}(i,j) = pvalue2;
                            pac_norm{5}(i,j) = pvalue3;
                            
                            
                            
                            if pvalue1<0.05 && pvalue2>0.05 && pvalue3 >0.05                          
                                phase_su{sigacc} = pfPha{1};
                                phase_unsu{sigacc} = pfPha{2};
                                sig_f{sigacc} = f1(i);                                
                                whatRcomputing2 = string(strcat(base,'_',other_ROIs,hemi));                                
                                sig_thiscon{sigacc} = whatRcomputing2;
                                MIs{sigacc} = ([M_norm{1}]);
                                MIuns{sigacc} =([M_norm{2}]);
                                sigacc= sigacc+1;
                                
                            else
                                insig_phase_su{insigacc} = pfPha{1};
                                insig_phase_unsu{insigacc} = pfPha{2};
                                insig_f{insigacc} = f1(i);
                                whatRcomputing2 = string(strcat(base,'_',other_ROIs,hemi));
                                insig_thiscon{insigacc} = whatRcomputing2;
                                insig_MIs{insigacc} = ([M_norm{1}]);
                                insig_MIuns{insigacc} =([M_norm{2}]);
                                insigacc= insigacc+1;
                                
                                                                
                            end

                            
                            
                        end
                    end
                    toc
                    
                    
                    path2save = string(strcat(cd,'/xPAC_surromax/',base,'-',other_ROIs,hemi,'/'));
                    if ~exist(path2save,'dir')
                        mkdir(path2save)
                    end
                    whatRcomputing = string(strcat(base,'_',other_ROIs,hemi,num2str(paircount)));
                    sprintf(whatRcomputing)
                    file2save1 = strcat(path2save,subjects,'_',whatRcomputing,'.mat');
                    save(string(file2save1),'pac_norm');
                    paircount = paircount+1;
                    

                    
                    
                    insig_MI = [insig_f;insig_thiscon;insig_MIs;insig_MIuns];
                    insig_phaamp = [insig_phase_su;insig_phase_unsu];
                    path2save2 = string(strcat(cd,'/xPAC_surromax/Insig_pairs/'));
                    if ~exist(path2save2,'dir')
                        mkdir(path2save2)
                    end
                    file2save2 = strcat(path2save2,other_ROIs,'_MI.mat');
                    file2save3 = strcat(path2save2,other_ROIs,'_phaamp.mat');
                    save(string(file2save2),'insig_MI')
                    save(string(file2save3),'insig_phaamp')
                    
                    
                    if sigacc >1
                        sig_MI = [sig_f;sig_thiscon;MIs;MIuns];
                        sig_phaamp = [phase_su;phase_unsu];
                        % save 2
                        path2save2 = string(strcat(cd,'/xPAC_surromax/Sig_pairs/'));
                        if ~exist(path2save2,'dir')
                            mkdir(path2save2)
                        end
                        file2save2 = strcat(path2save2,other_ROIs,'_MI.mat');
                        file2save3 = strcat(path2save2,other_ROIs,'_phaamp.mat');
                        save(string(file2save2),'sig_MI')
                        save(string(file2save3),'sig_phaamp')
                    end
                    
                    
                end
            end
        end
    end
    
end



sprintf('finished')