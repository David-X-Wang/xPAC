function Nochannel = get_channelno(ROIs,subject)
% Given ROIs and subjects NO. to find ROIs channel numbers of the subject
% ROIs is a 1*n string cell that each cell contains one ROI. 
%  Nochannel = get_channelno(ROIs,subject)
% subject is a 1*m string cell that each cell contains one. i.e. 'UT045'
% Nochannels  is a m*n string cell.
% Example: 
%   ROIs = {'AH left','AH right','PH left','PH right'}
%   subjects = {'UT047','UT053'}
%   Nochannels returns a 2*4 cell indicated ROIs channel number of each
%   subjet.
% this function works with info.mat which is converted from Brodmann_areas_new1.xlsx
load info.mat
for i = 1:length(ROIs)
    for j = 1:length(subject)
        index_y(i) = find(strcmp({info{1,:}},ROIs{i})==1);
        index_x(j) = find(strcmp({info{:,1}},subject{j})==1);
        Nochannel{j,i}= [info{index_x(j),index_y(i)}];
    end
end
end
