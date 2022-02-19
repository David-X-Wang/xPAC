
function Nochannel = whohasROI(subject,ROIs)

load info.mat
load ROIabbre.mat

index_y = zeros(1,length(ROIs));
index_x = zeros(1,length(subject));
%Nochannel = cell(length(ROIs),length(subject));
for i = 1:length(ROIs)
    for j = 1:length(subject)
        
        if isempty(find(strcmp({info{1,:}},ROIs{i})==1))
            Nochannel{i,j} = NaN;
            continue
        end
        index_y(i) = find(strcmp({info{1,:}},ROIs{i})==1);
        index_x(j) = find(strcmp({info{:,1}},subject{j})==1);
        Nochannel{i,j} = [info{index_x(j),index_y(i)}];
        if isempty(Nochannel{i,j})==0
            subj(i,j) = subject(j);

        end
    end
    
end
%subjhas = {ROIs' subj};
end
