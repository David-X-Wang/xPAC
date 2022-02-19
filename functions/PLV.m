function PLVnorm = PLV(Base_phase,Other_phase)



PLVnorm = zeros(1,size(Base_phase,1));


for i = 1:size(Base_phase,1)
    
    roipha = Other_phase(i,:);
    phase = Base_phase(i,:);
    
    
    zn = zeros(1,size(phase(:,21:end),2));
    for  j = 1:size(phase(:,21:end),2)
        zn(j) = exp(1i*(phase(20+j)-roipha(20+j)));
    end
    
     PLVnorm(i) = abs(mean(zn));
end
end