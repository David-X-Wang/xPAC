function M_raw = xPAC(Base_phase,Other_amp)

for i = 1:size(Base_phase,1)
    for  j = 1:size(Base_phase(:,51:end),2)
        zn(i,j) = Other_amp(i,50+j)*exp(1i*Base_phase(i,50+j));
    end
end
    M_raw = abs(mean(zn'));
end