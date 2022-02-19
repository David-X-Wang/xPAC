


subjects_total =  {'UT021','UT027','UT049','UT050','UT053','UT056','UT063','UT064','UT067','UT069',...
    'UT074','UT083','UT084','UT090','UT095','UT113','UT114','UT117','UT122','CC013','CC014','CC015',...
    'CC016','CC017','UT008','UT009','UT011','UT013','UT018','UT020','UT021','UT025','UT034',...
    'UT037','UT039','UT077','UT079','UT081','UT103','UT111'};


for i =1:length(subjects_total)
    
    file2load = sprintf('/users/xiaoliangw/EEG_DATA/%s.mat',subjects_total{i});
    eventstotal = ([events.recalled]);
    Noevents_all = length(eventstotal);
    Noevents_recalled = length(find(eventstotal==1));
    recallrate(i) = Noevents_recalled/Noevents_all;
    listtotal(i) = double(max([events.list]));
end

mean(recallrate)
std(recallrate)


