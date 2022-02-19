% coordinates visualization
clc 
clear
close all

addpath('/Volumes/Study/PAC/PACwithinSite/Visualization/Corrected/')
addpath('/Volumes/Study/')
addpath('/Volumes/Study/Matlab_tools/BrainNetViewer_20181219/')
% Left


%%
% Load node files from two hemisperes
% Load Edge Files from two hemispheres


% load AH_L_slow.node
% load AH_R_slow.node
% node = [AH_L_slow; AH_R_slow];
% save('AH_slow.node','node','-ascii','-tabs')
%
% load AH_L_slow.edge
% load AH_R_slow.edge
% edge = blkdiag(AH_L_slow, AH_R_slow);
% save('AH_slow.edge','edge','-ascii','-tabs')
%
% BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv','AH_slow.node','AH_slow.edge','allsig.mat');

%title('Significant Elec. Pairs of xPAC between Su. and Unsu. Encoding, AH_\theta(slow)  ROI_\gamma(slow)')
% title('Significant Elec. Pairs of xPAC during encoding, AH_\theta(slow)  ROI_\gamma(slow)')


%basetotal = {'AH','PH'}; % base region that theta phase comes from, chosen from AH or PH
bandtotal = {'slow','fast','slowfast','fastslow'};

%for b = 1:length(basetotal)
    for c = 1:length(bandtotal)
%        base = basetotal{b};
        band = bandtotal{c};
        
        whatRcomputing = strcat('_',band);
        
        load(strcat('L_',band,'.node'))
        load(strcat('R_',band,'.node'))
        
%        whatRcomputing = strcat(base,'_',band);
%         if exist(strcat(base,'_L_',band,'.node'))
%         load(strcat(base,'_L_',band,'.node'))
%         nodeL = eval(strcat(base,'_L_',band));
%         else
%             nodeL=[];
%         end
%         if exist(strcat(base,'_R_',band,'.node'))
%         load(strcat(base,'_R_',band,'.node'))
%         nodeR = eval(strcat(base,'_R_',band));
%         else
%             nodeR=[];
%         end
%         
%         node = [nodeL; nodeR];
%         save(strcat('/Volumes/Study/PAC/PLV/Visualization/Corrected/',...
%             base,'_',band,'.node'),'node','-ascii','-tabs')
        
        node = [eval(strcat('L_',band)); eval(strcat('R_',band))];
        
        save(strcat('/Volumes/Study/PAC/PACwithinSite/Visualization/Corrected/',...
            band,'.node'),'node','-ascii','-tabs')
        
        
%         load(strcat(base,'_L_',band,'.edge'))
%         load(strcat(base,'_R_',band,'.edge'))
%         edge = blkdiag(eval(strcat(base,'_L_',band)), eval(strcat(base,'_R_',band)));
%         save(strcat('/Volumes/Study/PAC/xPAC_ele_level/Visualization/Corrected/',...
%             base,'_',band,'.edge'),'edge','-ascii','-tabs')
        
        
        BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv',strcat(band,'.node'),'allsig.mat');
        
%        bb = strcat(base,'_\theta  ',band(1:4));
%         
%         if length(band)<5
%             roii = strcat('ROI_\gamma  ',band(1:4));
%         else
%             roii = strcat('ROI_\gamma  ',band(5:8));
%         end
        
       title(strcat('Significant Elec. Pairs of PAC during encoding..',band),'fontsize',16)
       saveas(gcf,strcat('/Volumes/Study/PAC/PACwithinSite/Visualization/Corrected/',band,'.png'))
        
    end
%end