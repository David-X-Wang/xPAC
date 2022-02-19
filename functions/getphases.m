% get phases for interests.
% [recallPhase, nonrecallPhase] =getphases(subjects,phaseROI,hemi,thetaOI)

function [recallPhase, nonrecallPhase] =getphases(subjects,phaseROI,hemi,thetaOI)






ROIs = string(strcat(phaseROI,'_',hemi));

ROI = ROIabbr(ROIs);



Nochannel =whohasROI(subjects,ROI);
acc =1;

for i = 1:length(Nochannel)
    chan_indx =  str2num(Nochannel{i});
    
    load(string(subjects{i}));
    
    index_recalled = find([events.recalled]==1);
    index_nonrecalled = find([events.recalled]==0);
    recalled_eeg = eeg(:,index_recalled,:);
    nonrecalled_eeg = eeg(:,index_nonrecalled,:);
    
    
    base_ind_total = find(ismember(channels, chan_indx));
    
    for base_ind = 1:length(chan_indx)
        ah{1} = reshape(recalled_eeg(base_ind_total(base_ind),:,:),[size(recalled_eeg,2),size(recalled_eeg,3)]);
        ah{2} = reshape(nonrecalled_eeg(base_ind_total(base_ind),:,:),[size(nonrecalled_eeg,2),size(nonrecalled_eeg,3)]);
        
        
        theta{1} = bandpass(ah{1}',thetaOI,250);
        theta{2} = bandpass(ah{2}',thetaOI,250);
        
        theta{1} = theta{1}';
        theta{2} = theta{2}';
        
%          z_theta{1} = zeros(size(theta{1}));
%          z_theta{2} = zeros(size(theta{2}));
         
        for j = 1:size(theta{1},1)
            z_theta{1}(j,:) = hilbert(theta{1}(j,:));

        end
        
        for j = 1:size(theta{2},1)
        z_theta{2}(j,:) = hilbert(theta{2}(j,:));
        end
        
        thetaphase{1} = angle(z_theta{1});
        thetaphase{2} = angle(z_theta{2});
        gammaamp{1} = abs(z_theta{1});
        gammaamp{2} = abs(z_theta{2});
        
        recallPhase{acc}=thetaphase{1};
        nonrecallPhase{acc}=thetaphase{2};
        
        

        acc = acc+1;
        clear z_theta
    end
           
end

end

% ROI_ind_total=find(ismember(channels,ROI_indx));
% 
% base_ind = 1;
% ROI_ind =1;
% 
% 
% 
% para{1} = reshape(recalled_eeg(ROI_ind_total(ROI_ind),:,:),[size(recalled_eeg,2),size(recalled_eeg,3)]);
% para{2} = reshape(nonrecalled_eeg(ROI_ind_total(ROI_ind),:,:),[size(nonrecalled_eeg,2),size(nonrecalled_eeg,3)]);



