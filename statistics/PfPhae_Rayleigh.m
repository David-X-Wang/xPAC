%Prefered phase test
% Rayleigh test for each prefered phase against zero(uniform distribution)
% and Watson-Williams test for each two pair against each other (Left-right,
% and Anterior-Posterior)
clear
close all
clc


addpath('/Volumes/Study/PAC/statistics/CircStat');

path1 = ('/Volumes/Study/PAC/f-f_inspection/xPAC/xPAC_Surro/Sig_pairs/');

addpath(path);

ROIs = {'EC','paraHG'};%'PC','LPF','SM','LMT'};

con_acc= 1;
for rind = 1:length(ROIs)
    
    
    filename1 = (strcat(ROIs(rind),'_MI.mat'));
    filename2 = (strcat(ROIs(rind),'_phaamp.mat'));
    
    load(string(strcat(path1,filename1)))
    load(string(strcat(path1,filename2)))
    
    conn = unique([sig_MI{2,:}]);
    
    figure('position',[0 0 600 200])
    for i = 1:length(conn)
        %connections = unique([totalMI{2,:}]);
        
        phaind1 = find(ismember([sig_MI{2,:}],conn(i)));
        %phaind2 = find(ismember([totalMI{2,:}],conn(i)));
        A = sig_phaamp(1,phaind1);
        %B = totalphaamp(1,phaind2);
        
        for j = 1:length(A)
            [p{i}(j),~] =  circ_rtest(A{j});
        end
        groupind{i} = j*ones(size(A));
        
        groupname{i} = conn{i};
        
    end
    
    boxplot(cell2mat(p),cell2mat(groupind))
    ylabel('P-Value')
    set(gca, 'XTickLabel', groupname,'fontsize',14);
    title('Rayleigh Test of Phases against uniform dist.')
    clear p groupind groupname
    saveas(gca,string(strcat(path1,'Rayleightest_',ROIs(rind),'.png')))
end
%% hemi ttest







