function b = ROIabbr(a)

load ROIabbre
ind =zeros(1,length(a));
for i = 1:length(a)
ind(i) = find(strcmp({ROIabbre{1,:}},a{i})==1);
end
b = {ROIabbre{2,ind}};

end
