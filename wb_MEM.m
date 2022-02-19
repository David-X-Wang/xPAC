% wide-band xPAC connectivity analysis
% this script computes the z-score of certian xPAC/PLV over shuffled
% distribution (considered as pdf of the connectivity), compared with sample
% distribution.


% data matrix cell form: [slow-slow,slow-fast;fast-slow,fast-fast]

clc
clear all
close all


addpath('/Volumes/Study/PAC/statistics/');
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;            % 1-tailed t-distribution
% T-Statistic Given Probability alpha & Degrees-Of-Freedom v
t_inv=@(a,v) fzero(@(tval) (max(a,(1-a))-tdist1T(tval,v)),5);
alpha_to_Z = @(a) -sqrt(2) * erfcinv(a*2);%get Z value given the alpha

% path2load = ('/Volumes/Study/PAC/wide-band/wb/xPAC');
% path2save = ('/Volumes/Study/PAC/wide-band/wb/memodel/behav_network/');


path2load = ('/Volumes/Study/PAC/wide-band/wb/xPAC');
path2save = ('/Volumes/Study/PAC/wide-band/wb/xPAC/');


baset= {'AH','PH'};
hemit = {'_L','_R'};
ROIs = {'EC','paraHG','LPF','LMT','PC','BTL','SM'};


labelL = strcat([baset,ROIs],hemit(1));
labelR = strcat([baset,ROIs],hemit(2));
Labels = [labelL',labelR'];


cal = 'xPAC';
if strcmp(cal,'xPAC')
    band = {'slow','fastslow','slowfast','fast'};
else
    band = {'slow','fast'};
end

acc= 1;


Z_nonbehav= cell(1,2);
Z_nonbehav{1} = zeros(9,9);
Z_nonbehav{2} = zeros(9,9);

Z_behav= cell(1,2);
Z_behav{1} = zeros(9,9);
Z_behav{2} = zeros(9,9);

% myratio = cell(1,length(band));

for bandind = 1:length(band)
    
    FixFX= cell(1,2);
    Intercepts = cell(1,2);
    Randeffect = cell(1,2);
    for hind = 1:length(hemit)
        
        intercept = zeros(2,7);
        fixedeffect =  zeros(2,7);
        fixedeffect_T = zeros(2,7);
        rand_effect = zeros(2,7);
        
        for baseind = 1:length(baset)
            
            base = baset(baseind);
            
            for Rind = 1:length(ROIs)
                
                
                
                hemi = hemit(hind);
                ROI = ROIs(Rind);
                
                file2load1 = sprintf(strcat(path2load,'/%s.mat'), string(strcat( base,'_',ROI,hemi)));
                
                load(file2load1)
                
               % file2load2 = sprintf(strcat(path2load,'/%s_stats.mat'), string(strcat( base,'_',ROI,hemi)));
                
                data_su = MI_total(1,:);
                data_unsu = MI_total(2,:);
                %data_sunsu = stats_total(1,:);
                
                MIs_su = cellfun(@(v)v(bandind),data_su);
                MIs_unsu = cellfun(@(v)v(bandind),data_unsu);
                NoEventsu= cellfun(@length,MIs_su);
%                 MIs_su = cellfun(@nanmean,MIs_su, 'UniformOutput', false);
%                 MIs_unsu = cellfun(@nanmean,MIs_unsu, 'UniformOutput', false);
                
                N=cellfun(@length,MIs_su);
                
                %                 for i = 1:length(MIs_su)
                %                     MIs_unsu{i} = MIs_unsu{i}(randi(length(MIs_unsu{i}),1,N(i)));
                %                 end
                
                
               % subjNo = 1:length(unique(NoEventsu));
                NoEventunsu= cellfun(@length,MIs_unsu);
                
                Noeles= cellfun(@length,MIs_su);
                Noeles_unsu = cellfun(@length,MIs_unsu);
                
                
%                 elepairsu = cell(1,length(Noeles));
%                 elepairunsu = cell(1,length(Noeles));
                
                
                MIs_su=[MIs_su{:}];
                MIs_unsu = [MIs_unsu{:}];
                
                MIs=[ MIs_unsu MIs_su];
                 condi = cell(1,length(MIs));
                condi(1:length(MIs_unsu))={'unsu'};
                condi(1+length(MIs_unsu):end)={'su'};
                
                
%                 condi = zeros(1,length(MIs));
%                 condi(1:length(MIs_unsu))=0;
%                 condi(1+length(MIs_unsu):end)=1;
                
                subjNeventsu = cell(1,length(NoEventsu));
                subjNeventunsu = cell(1,length(NoEventunsu));
                
                for i=1:length(NoEventsu)
                    subjNeventsu{i} = cellstr(num2str(NoEventsu(i)*ones(Noeles(i),1)))';
                    subjNeventunsu{i} = cellstr(num2str(NoEventsu(i)*ones(Noeles_unsu(i),1)))';
                end
                
                subjNeventsu = [subjNeventsu{:}];
                subjNeventunsu = [subjNeventunsu{:}];
                subjs = [subjNeventunsu subjNeventsu];
                % tbl=table(MIs_su',MIs_unsu',subjNevents','variablename',{'MIsu','MIunsu','subjNevents'});
                tbl=table(MIs',condi',subjs','variablename',{'MI','condi','subj'});
                lme = fitlme(tbl,'MI~condi+(1|subj)');
                
                %  lme = fitlme(tbl,'MIsu~MIunsu+(1|subjNevents)');
                
                
                
                intercept(baseind,Rind) = double(lme.Coefficients(1,4));
                fixedeffect(baseind,Rind) = double(lme.Coefficients(2,6));
                fixedeffect_T(baseind,Rind) = double(lme.Coefficients(2,4));
                rand_effect(baseind,Rind) = std(randomEffects(lme));
                
                
                
                
                
                %Z{baseind,Rind}
                
                %  Z{baseind,Rind} =  z_behav;
                
            end
            
        end
        
        
        Intercepts{hind}(1:2,3:9) = intercept;
        Intercepts{hind}(3:9,1:2) = intercept';

        

%        p = 0.5 * erfc(fixedeffect_T ./ sqrt(2));
%       FixFX{hind}(1:2,3:9)= -alpha_to_Z(p);
%       FixFX{hind}(3:9,1:2) = -alpha_to_Z(p)';
        FixFX{hind}(1:2,3:9) = fixedeffect_T;
        FixFX{hind}(3:9,1:2) = fixedeffect_T';
        
   
       
     
       
        
        % BehavFXcorrected = ;
    end
    
   
     
     
    Intercepts = blkdiag(Intercepts{1},Intercepts{2});
    Randeffect = blkdiag(Randeffect{1},Randeffect{2});
    
     FixFX = blkdiag(FixFX{1},FixFX{2});
    p = 0.5 * erfc(FixFX ./ sqrt(2));
    p(p==0.5)=0; 
    [~,~,~,BehavFXp] = fdr_bh(p); 

    BehavFXcorrected = -alpha_to_Z(BehavFXp);
    BehavFXcorrected(BehavFXcorrected==Inf)=0;
    
%     BehavFXcorrected = FixFX;
%     BehavFXcorrected(BehavFXcorrected==Inf)=0;
   
    
    BehavFXcorrected(BehavFXcorrected>3)=2.999;
    BehavFXcorrected(BehavFXcorrected<-3)=-2.999;
    BehavFXcorrected(BehavFXcorrected~=0)=BehavFXcorrected(BehavFXcorrected~=0)+3;
  
    BehavFXcorrected(BehavFXcorrected<3+1.64&BehavFXcorrected>3-1.64)=0;

    BehavFXcorrected(1,1)=6;

    BehavFXcorrected=BehavFXcorrected*100;
    BehavFXcorrected(2,2)=1; 
    
    
  %       path2save = '/Volumes/Study/PAC/wide-band/wb/memcorrection/';
    if ~exist(path2save,'dir')
        mkdir(path2save)
    end
    
    file2save1 = sprintf(strcat(path2save,'/MEM3_%s.csv'), band{bandind});
    csvwrite(file2save1,BehavFXcorrected)
    
%     file2save2 = sprintf(strcat(path2save,'/intercept_%s.csv'), band{bandind});
%     csvwrite(file2save2,Intercepts)
%     
%     file2save3 = sprintf(strcat(path2save,'/randeffect_%s.csv'), band{bandind});
%     csvwrite(file2save3,Randeffect)
end


