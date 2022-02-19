%% fdr correction

clear 
close 
clc

addpath('/Volumes/Study/PAC/xPAC/Pvalue_results/Ca')
hemi = {'_L','_R'}; % hemisphere chosen from _L and _R
base = {'AH','PH'}; % base region that theta phase comes from, chosen from AH or PH
band = {'slow','fast','slowfast','fastslow'};


for a = 1:length(hemi)
    for b = 1:length(base)
        for c = 1:length(band)
            
            
            whatRcomputing = string(strcat(base(a),hemi(b),'_',band(c)));
            sprintf(whatRcomputing)
            file2load = strcat('/Volumes/Study/PAC/xPAC/Pvalue_results/Ca/',whatRcomputing,'_pvalue.mat');
            
            load(file2load)

            for i =1:1
                for j=1:8
                    fdr{i}{j} = mafdr(PA{i}{j});
                end
            end
            
            PA = fdr';
            clear fdr
            file2save = strcat('/Volumes/Study/PAC/xPAC/Pvalue_results/Corrected/',whatRcomputing,'_pvalue.mat');
            save(string(file2save),'PA');
        end
    end
end


