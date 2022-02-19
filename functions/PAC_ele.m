% Cross-Region Phase Amplitude Coupling at electrodes level
% This script computes the xPAC of base region to other regions
% This script extracts ROIs of Subjects and will return
% # channel of ROIs of given subjects. the Nochannel is
% a #channel by # subjects matrix.
% David Wang 03/20/2019

function xMI = PAC_ele(subject,Nochannel,ROIs,band,Fs)

condition_total = {'successful','unsuccessful','retrieval'};
filename = strcat(subject,'.mat');
load(string(filename))

index_recalled = find([events.recalled]==1);
index_nonrecalled = find([events.recalled]==0);
recalled_eeg = eeg(:,index_recalled,:);
nonrecalled_eeg = eeg(:,index_nonrecalled,:);

filename = strcat(subject,'.mat');
load(string(filename))

Other_amp = cell(1,length(ROIs));

xMI = cell(length(condition_total),length(ROIs)-1);

for Cond = 1:length(condition_total)
    
    for Rind = 1:length(ROIs)
        
        % skip when subject has no related ROI
        if isempty(Nochannel{Rind}) ==1
            
            if Rind == 1
                Base_phase =[];
                continue
            else
                Other_amp{Rind} = [];
                continue
            end
        end
        
        tag_no = str2num(Nochannel{Rind});
        index = find(ismember(channels, tag_no));
               
        if isempty(index)
            continue
        end
                
        % forming EEG matrices 
        successful = recalled_eeg(index,:,:);  
        unsuccessful = nonrecalled_eeg(index,:,:); 
        retrieval = eeg_retrieval(index,:,:);
        
        condition = eval(condition_total{Cond});
        
        % preprocessing EEG of certain regions
        amp_gamma = cell(1,size(condition,1));
        phase_theta = cell(1,size(condition,1));   
        
        
        for i = 1:size(condition,1)
            x = reshape(condition(i,:,:),[size(condition,2),size(condition,3)]);
            [amp_gamma{i},phase_theta{i}]=PACprep(x,band,Fs);
        end
        
        Base_phase{Rind} = phase_theta;
        Other_amp{Rind} = amp_gamma;
        
%         if Rind == 1  % AH_Based   
%             Base_phase = phase_theta;
%         else          % if Rind == 2 % PH_Based
%             Other_amp{Rind-1} = amp_gamma;
%         end
          
    end
    
    % skip when subject has no Based region
    if isempty(Base_phase)
        continue
    end
    
    for Rind2 = 1:length(Other_amp)
        % skip when subject has no certain gamma region
        if isempty(Other_amp{Rind2}) ==1
            xMI{Cond,Rind2} = [];
            continue
        end     
        
        % cross region PAC
        acc = 1;
        M_store = cell(1,length(Base_phase)*length(Other_amp{Rind2}));
      
        for a = 1:length(Base_phase{Rind2})
            for b = 1:length(Other_amp{Rind2})
                M_raw = PAC(Base_phase{Rind2}{a},Other_amp{Rind2}{b});
                M_store{acc} = M_raw;
                acc = acc+1;
            end
            
        end       
        %M_raw = abs(mean(zn'));
        xMI{Cond,Rind2} = M_store;       
    end   
    clear condition
   % clear M_store
end
end

% subjMI{subjInd}=xMI;

% end


% file2save = strcat(cd,'/',base,hemi,'_',band,'.mat');
% %figtosave = strcat(path2save,'/',strcat(subject{:}),'_',strcat(condition_total(cond),'_',strcat(ROIa),'_low.fig'));
% save(string(file2save),'subjMI')


