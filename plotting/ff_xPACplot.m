% Fig3a f-f xPAC plots.
% David Wang
% 11/22/2019
%
clear;
close all;
clc;


df = 1;
df2 = 2;
f1 = 2:9;
f2 = 30:2:100;
thetalength = length(f1);
gammalength = length(f2);

hemit = {'_L','_R'};
baset = {'AH','PH'};
%ROIs = {'LPF','LMT','PC','BTL','SM'};

ROIs = {'EC','paraHG'};%,'BTL','SM','PC'};

    acc= 1;
figure('position',[0 0 1000 length(ROIs)*250])
for rind = 1:length(ROIs)
    ROI = ROIs{rind};

    for i = 1:length(hemit)
        hemi = hemit{i};
        for j = 1:length(baset)
            base = baset{j};                                                            
            file2load = strcat(base,'_',ROI,hemi);
                      
            load(sprintf('/Volumes/Study/PAC/f-f_inspection/ff/xPAC/%s.mat',file2load));
            
            
            conn_No=min(200,length(MI_total)); %number of files in the path
            Zsu = cell(1,conn_No);
            Zunsu = cell(1,conn_No);
            for cind= 1:conn_No
                Su =  cellfun(@mean, MI_total{1,cind});
                Unsu = cellfun(@mean, MI_total{2,cind});
                
                Zsu{cind} = Su;%(Su-mean(Su))./std(Su);
                Zunsu{cind} = Unsu;%(Unsu-mean(Unsu))./std(Unsu);
            end
            
            avgZsu = nanmean(cat(3,Zsu{:}),3);
            avgZunsu = nanmean(cat(3,Zunsu{:}),3);
                                 
            
            subplot(length(ROIs),4,acc)            
            imgsmooth(avgZunsu);
            set(gca,'YDir','normal')    
            
            
        
            if rind == length(ROIs)
                xlabel('Theta (Hz)');
            else
                ax2 = gca;
                set(ax2,'XTickLabel','');
            end
            
            if acc == 1
                ylabel('Gamma (Hz)');
            else
                ax2 = gca;
                set(ax2,'YTickLabel','')
            end
            n=50;
%             colormap(jet);
%             cmap = flipud(cbrewer('seq','OrRd',n)); % blues at bottom
%             colormap(cmap);
    
           % caxis([0.94 1.05])
            colorbar off
            if strcmp(ROI,'paraHG')
                roi = 'PHG';
            end
            titlename = string(strcat(base,'-',ROI,hemi));
            title(titlename)
            acc= acc+1;
            set(gca,'fontsize',12)
        end
    end
end

h = colorbar('position',[0.92 0.41 0.015 0.3]);

%print('-painters','-depsc','/Volumes/Study/PAC/f-f_inspection/ff_plotunsu.eps')

