
addpath('/Volumes/Study/PAC/plotting_tools');

clear all
close all
clc

hemit = {'_L','_R'};
baset = {'AH','PH'};
ROIs = {'EC','paraHG'};


acc=1;
for i = 1:length(hemit)
            MI = cell(1,2);
        band =cell (1,2);
    for j = 1:length(baset)
        
        acc=1;

        slowT = cell(1,2);
        fastT = cell(1,2);
        for rind = 1:length(ROIs)
            
            ROI = ROIs{rind};
            
            
            hemi = hemit{i};
            
            
            
            base = baset{j};
            file2load = strcat(base,'-',ROI,hemi);
            
            load(sprintf('/Volumes/Study/PAC/f-f_inspection/ff_mem/%s.mat',file2load));
           
            slow = Z(1:5,:);
            fast = Z(5:end,:);
            
            
            slowT{rind} = slow(:)';
            fastT{rind} = fast(:)';
            acc= acc+1;
        end
        
        slowT=[slowT{:}];
        fastT=[fastT{:}];
        MI{j} = [slowT fastT];
        band{j} = cell(1,length(MI{j}));
        band{j}(1:length(slowT)) = {'slow'};
        band{j}(1+length(slowT):end) = {'fast'};
    end
    base = cell(1,length(MI{1})+length(MI{2}));
    base(1:length(MI{1})) = {'AH'};
    base(1+length(MI{1}):end) = {'PH'};
    MI = [MI{:}];
    band = [band{:}];
    
    [~,~,stats] = anovan(MI,{base,band},'varnames',{'Base','Band'},'display',0);
    [c,mm,h,nms] = multcompare(stats,'Dimension',[1,2],'display','off');
    figure
    bar(mm(:,1))
    hold on
    errorbar(mm(:,1)',mm(:,2)')
   % ylim([-3 6])
    
    
    
end



