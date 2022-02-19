% Cross-Region Phase Amplitude Coupling function 
% This script computes the xPAC of base region to other regions 
% David Wang 03/10/2019

function xMI = xPACMI(subject,Nochannel)
Fs = 250;
band = 'fast';



if strcmp(band,'slow')
    gamma = [30 70]*2/Fs;
    theta = [2 5]*2/Fs;
else
    gamma = [70 110]*2/Fs;
    theta = [5 9]*2/Fs;
end

%%
condition_total = {'successful','unsuccessful','retrieval'};
% %subjects={'UT021'};
% filename = strcat(subject,'.mat');
% load(string(filename))
% gamma = [30 70]*2/Fs;
% theta = [2 5]*2/Fs;
or_theta =6;
or_gamma =12;
[b_gamma,a_gamma] = butter(or_gamma,gamma);
[b_theta,a_theta] = butter(or_theta,theta);

%%
% subjMI = cell(1,length(subject));
% for subjInd = 1:length(subjects)
    
   % subject = subjects(subjInd);
    filename = strcat(subject,'.mat');
    load(string(filename))
    
    Other_amp = cell(1,length(ROIs)-1);
    xMI = cell(3,length(ROIs)-1);
    
    for Cond = 1:3
        
        for Rind = 1:length(ROIs)
            
            if isempty(Nochannel{Rind}) ==1
                if Rind == 1
                    Base_phase =[];
                    continue
                else
                    Other_amp{Rind-1} = [];
                    continue
                end
            end
            
            tag_no = str2num(Nochannel{Rind});
            tag_no = tag_no(1);
            
            
            recall = zeros(1,size(eeg,2));
            
            for i = 1:size(eeg,2)
                recall(i) = events(i).recalled;
                %  word{i} = events(i).item;
            end
            
            index_recalled = find(recall==1);
            %index_recalled = index_recalled(index_recalled<=size(eeg,2));
            recalled = eeg(:,index_recalled,:);
            % find channel No. of ROI
%             index = 0;
            
%             if find(channels==tag_no(i))
                
            index = find(channels==tag_no);
            if isempty(index)
                continue
            end
            
            index_nonrecalled = find(recall==0);
            % index_nonrecalled = index_nonrecalled(index_nonrecalled<=size(eeg,2));
            nonrecalled = eeg(:,index_nonrecalled,:);
            % form three new matrix
            successful = recalled(index,:,:);  % recalled encoding ROI
            % retrieval = eeg_retrieval(index,:,:); % retrival ROI
            unsuccessful = nonrecalled(index,:,:); % non_recalled encoding ROI
            retrieval = eeg_retrieval(index,:,:);
            
            
            condition = eval(condition_total{Cond});
            
            x = reshape(condition,[size(condition,2),size(condition,3)]);
            
            x_gamma = zeros(size(x));
            x_theta = zeros(size(x));
            z_gamma = zeros(size(x));
            z_theta = zeros(size(x));
            phase_gamma = zeros(size(x));
            phase_theta = zeros(size(x));
            
            for i = 1:size(x,1)
                x_gamma(i,:) = filter(b_gamma,a_gamma,x(i,:));
                x_theta(i,:) = filter(b_theta,a_theta,x(i,:));
            end
            
            % gamma band operation
            
            for i = 1:size(x,1)
                z_gamma(i,:) = hilbert(x_gamma(i,:));
                z_theta(i,:) = hilbert(x_theta(i,:));
                phase_theta(i,:) = angle(z_theta(i,:));   % phase of theta
            end
            amp_gamma = abs(z_gamma);       % amplitdue of gamma
            
            if Rind == 1  % AH_Based
                % if Rind == 2 % PH_Based
                Base_phase = phase_theta;
            else
                Other_amp{Rind-1} = amp_gamma;
            end
            
        end
        
        
        for Rind2 = 1:length(Other_amp)
            if isempty(Other_amp{Rind2}) ==1
                xMI{Cond,Rind2} = [];
                continue
            end
            
            for i = 1:size(Base_phase,1)
                for  j = 1:size(Base_phase(:,51:end),2)
                    zn(i,j) = Other_amp{Rind2}(i,50+j)...
                        *exp(1i*Base_phase(i,50+j));
                end
            end
            M_raw = abs(mean(zn'));
            xMI{Cond,Rind2} = M_raw;
            
        end
%         clear condition
%         clear x
%         clear zn
    end
   % subjMI{subjInd}=xMI;
end
% end


% file2save = strcat(cd,'/',base,hemi,'_',band,'.mat');
% %figtosave = strcat(path2save,'/',strcat(subject{:}),'_',strcat(condition_total(cond),'_',strcat(ROIa),'_low.fig'));
% save(string(file2save),'subjMI') 


