% test analytic signal for circuarity by generalized likelyhood ratior
% test(GLRT). 
% function for script sigCircu.m
% 09/19/2019

%Z = randn(5,100)+1i*randn(5,100);

function rho = circoeff(Z)

rho = zeros(1,size(Z,1));
for i = 1:size(Z,1)
    z = Z(i,:);
    var = mean(abs(z).^2);
    q = mean(z.^2)/var;
    rho(i)= abs(q);
end
end