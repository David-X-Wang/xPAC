function [M_raw,Z,pfPha] = xPAC(Base_phase,Other_amp)

pfPha = zeros(1,size(Base_phase,1));
M_raw = zeros(1,size(Base_phase,1));


for i = 1:size(Base_phase,1)
    
    amplitude = Other_amp(i,:);
    phase = Base_phase(i,:);
    
    zn = zeros(1,size(phase(:,21:end),2));
    for  j = 1:size(phase(:,21:end),2)
        zn(j) = amplitude(20+j)*exp(1i*phase(20+j));
    end
    pfPha(i) = angle(mean(zn));    
    M_raw(i) = abs(mean(zn));
    Z(i,:) = zn;   
end   
end