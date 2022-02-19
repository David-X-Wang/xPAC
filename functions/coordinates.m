% coordinates extraction
% This script extracts the xyz coordinates of the significant electrode
% pairs of xPAC.
% V0.1
% David Wang
% lastest update: 04/25/2019
clc;
clear;
close;



subjects =  {'UT021','UT027','UT049' 'UT050','UT053',...
    'UT056','UT063','UT064', 'UT067','UT069',...
    'UT074','UT083','UT084','UT090','UT095',...
    'UT113','UT114','UT117','UT122'};

load info.mat
load ROIabbre.mat


hemit = {'L','R'}; % hemisphere chosen from _L and _R
bandt = {'slow','fast','slowfast','fastslow'};
other_ROIs = {'EC','paraHG','LPF','LMT','PC','BTL','SM'};



for a = 1:length(hemit)
    
    for c = 1:length(bandt)
        
        hemi = hemit{a}; % hemisphere chosen from _L and _R
        % base region that theta phase comes from, chosen from AH or PH
        band = bandt{c};
        
        ROIs  = ['AH','PH', other_ROIs];
        
        
        
        ROIs =ROIabbr(strcat(ROIs,'_',hemi)); % left hemisphere
        
        [subjhas, Nochannel] = whohasROI(subjects,ROIs);
        
        nsub = length(subjects);
        
        whatRcomputing = strcat('_',hemi,'_',band);
        sprintf(whatRcomputing)
        
        file2load = strcat(cd,'/Pvalue/Corrected/',whatRcomputing,'_pvalue.mat');
        
        load(string(file2load));
        
        
        
        for j = 1: size(Nochannel,1) % ROI index
            acc = 1;
            
            for i = 1:nsub    % subject index
                %      base_chan = str2num(Nochannel{1,i});
                
                roi_chan = str2num(Nochannel{j,i});
                
                len(j,i) = length([roi_chan]);%*length([roi_chan]);
                %leng = leng+len;
                %  for k = 1:length(base_chan)
                
                %             if isempty (base_chan(k))
                %                 continue
                %             end
                
                for l = 1:length(roi_chan)
                    if isempty (roi_chan(l))
                        continue
                    end
                    channo{j,acc} = [roi_chan(l)];
                    
                    acc = acc+1;
                end
                % end
                
            end
        end
        
        
        
        sig_chann=cell(size(channo));
        sig_chann_p = ones(size(channo));
        
        for j = 1: size(Nochannel,1) % ROI index
            for i = 1:size(Nochannel,2)   % subject index
                
                left = sum(len(j,1:i-1))+1;
                right = sum(len(j,1:i));
                
                if right> length(PA{1}{j})
                    continue
                end
                %   sig_ind= find(PA{1}{j-1}(left:right)<0.05);
                sig_ind= left:right;
                if isempty(sig_ind)
                    continue
                end
                %         sig_chan_p{j-1} = PA{1}{j-1}(left:right)<0.05;
                for k = 1:length(sig_ind)
                    sig_chann{j,sig_ind(k)} = channo{j,sig_ind(k)};
                    sig_chann_p(j,sig_ind(k)) = PA{1}{j}(sig_ind(k));
                end
                
            end
        end
        
        
        
        
%         for i = 1:size(sig_chann,1)
%             for j = 1:size(sig_chann,2)
%                 if sig_chann_p(i,j)>0.05
%                     sig_chann{i,j}=[];
%                 end
%             end
%         end
        sig_chann_p(sig_chann_p>0.05)=1;
        %%
        sigsub = cell(1,size(len,2));
        sigsub_p = cell(1,size(len,2));
        for i = 1:size(len,2)
            for j = 1:size(len,1)
                left = sum(len(j,1:i-1))+1;
                right = sum(len(j,1:i));
                zz = sig_chann(j,left:right);
                xx = sig_chann_p(j,left:right);
                sigsub{i}{j,:} = zz(~cellfun(@isempty,zz));
                sigsub_p{i}{j,:} = xx;%(xx~=0);
            end
        end
        
        
        %% extract coordinates from the talStract file
        coordi = sigsub;
        total_coordi = cell(1,length(sigsub));
        total_p = cell(1,length(sigsub));
        
        acc2 =1;
        
        for sub = 1:length(sigsub)
            if cellfun(@isempty, sigsub{sub})
                continue
            end
            
            filename = strcat(subjects{sub},'.mat');
            load(string(filename))
            
            filename2 = strcat('/Volumes/Study/EEG_DATA/',subjects{sub},'_talLocs_database_monopol.mat');
            if isfile(filename2)
                load(string(filename2))
            else
                continue
            end
            
            
            if exist('talStruct') == 0
                continue
            end
            
            coor = [talStruct.x; talStruct.y; talStruct.z ]';
            
            
            
            acc = 1;
            for roi = 1:length(sigsub{sub})
                
                if isempty(sigsub{sub}{roi})
                    continue
                end
                
                for i = 1:length(sigsub{sub}{roi})
                    ele1 = sigsub{sub}{roi}{i};
                    ele = find(ismember(channels,ele1))';
                    if isempty(ele)
                        continue
                    end
                    %             if length(ele)<2
                    %                 continue
                    %             end
                    
                    coordi{sub}{roi}{i} = coor(ele,:)';
                    total_p{sub} = [sigsub_p{sub}{:}];
                    total_coordi{sub}{acc,:} = coor(ele,:)';
                    node_gen(acc2,:) = [coor(ele,:), 2];
                    
                    acc = acc+1;
                    acc2 = acc2+1;
                end
            end
            clear talStruct
        end
        %%

        node_value = total_p(~cellfun(@isempty,total_p));
        node_value = [node_value{:}];
        node_value(isnan(node_value))=0.001;
        
        node_final = [node_gen node_value(1:length(node_gen))' ones(size(node_gen,1),1)];
        

        temp=node_final(:,5);
        node_final(:,5)=node_final(:,4);
        node_final(:,4)=temp;
        
        node_final(node_final(:,4)>0.05,5)=1;
        
        %edge =zeros(size(node_gen,1),size(node_gen,1));
        clear node_gen
        
        
        % edge_value = 2-edge_value;
        
        % for i = 1:acc2-1
        %     edge(i*2-1,i*2)=node_value(i);
        %     edge(i*2,i*2-1)=node_value(i);
        % end
        
        %% file saving
        
        whatRcomputing = strcat(hemi,'_',band);
        sprintf(whatRcomputing)
        file2save1 = strcat('/Volumes/Study/PAC/PACwithinSite/Visualization/Corrected/',whatRcomputing,'.node');
        %file2save2 = strcat('/Volumes/Study/PAC/xPAC_ele_level/Visualization/Corrected/',whatRcomputing,'.edge');
        save(string(file2save1),'node_final','-ascii','-tabs');
        %save(string(file2save2),'edge','-ascii','-tabs');
        
        
    end
end

