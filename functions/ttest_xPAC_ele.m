% 2-sample t-test of xPAC computation results
% [P1, P2 ,P3] = ttest_xPAC_ele(xMI_total)
%
%
% David Wang 03/20/2019


function [P, M, STD] = ttest_xPAC_ele(xMI_total)

% two-sample t-test
for subjInd = 1:length(xMI_total)
    P1{subjInd} = cell(1,8);
    P2{subjInd} = cell(1,8);
    P3{subjInd} = cell(1,8);
    M_recall{subjInd} = cell(1,8);
    M_nonrecall{subjInd} = cell(1,8);
    M_retrieval{subjInd} = cell(1,8);
    STD_recall{subjInd} = cell(1,8);
    STD_nonrecall{subjInd} = cell(1,8);
    STD_retrieval{subjInd} = cell(1,8);
    
    if isempty(xMI_total{subjInd})
        continue
    end
    
    
    for Rind = 1:length(xMI_total{subjInd})
        % successful encoding VS unsuccessful encoding
        if isempty(xMI_total{subjInd}{1,Rind})
            % P1{subjInd}{eind,Rind} = [];
            continue
        end
        
        
        for eind = 1:length(xMI_total{subjInd}{1,Rind})
            
            
            if isempty(xMI_total{subjInd}{1,Rind}{eind})
                P1{subjInd}{eind,Rind} = [];
                continue
            end
            
%             M_recall{subjInd}{eind,Rind} = nanmean(xMI_total{subjInd}{1,Rind}{eind});
%             M_nonrecall{subjInd}{eind,Rind} = nanmean(xMI_total{subjInd}{2,Rind}{eind});
%             M_retrieval{subjInd}{eind,Rind} = nanmean(xMI_total{subjInd}{3,Rind}{eind});
%             STD_recall{subjInd}{eind,Rind} = nanstd(xMI_total{subjInd}{1,Rind}{eind});
%             STD_nonrecall{subjInd}{eind,Rind} = nanstd(xMI_total{subjInd}{2,Rind}{eind});
%             STD_retrieval{subjInd}{eind,Rind} = nanstd(xMI_total{subjInd}{3,Rind}{eind});
            [~,P1{subjInd}{eind,Rind}] = ttest2(xMI_total{subjInd}{1,Rind}{eind},...
                xMI_total{subjInd}{2,Rind}{eind},'tail','both');
            
            
            %             if P1{subjInd}{eind,Rind} < 0.1
            %                 %                 su1{Rind}(acc) = mean(xMI_total{subjInd}{1,Rind}{eind});
            %                 %                 unsu1{Rind}(acc) = mean(xMI_total{subjInd}{2,Rind}{eind});
            %                 %             else
            %                 su2{Rind}(eind) = mean(xMI_total{subjInd}{1,Rind}{eind});
            %                 unsu2{Rind}(eind) = mean(xMI_total{subjInd}{2,Rind}{eind});
            %                 acc = acc+1;
            %             end
            
            
%             [~,P2{subjInd}{eind,Rind}] = ttest2(xMI_total{subjInd}{1,Rind}{eind},...
%                 xMI_total{subjInd}{3,Rind}{eind},'tail','both');
%             
%             
%             [~,P3{subjInd}{eind,Rind}] = ttest2(xMI_total{subjInd}{2,Rind}{eind},...
%                 xMI_total{subjInd}{3,Rind}{eind},'tail','both');
            
            
        end
    end
end
P = {P1;P2;P3} ;
M = {M_recall;M_nonrecall;M_retrieval};
STD = {STD_recall;STD_nonrecall;STD_retrieval};
end

% % Z-scorring
% for
%   zRec = P1;
%   zRec = -norminv(zRec);
% end
%
%   Z= [P1];

