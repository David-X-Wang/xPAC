% Preprocessing of signal to cross-region PAC. 
%
% this function filters input signal to gamma and theta band, 
% then returns the retures amplitude(envolope) of gamma band
% and phase of theta band by hilbert transform. 
% 
%   [amp_gamma,phase_theta]=PACprep(x,band,Fs);
%   x is input signal(EEG trials)  
%   Fs is the sampling frequency of the EEG data
%   band is chosen from 'slow' 2-5Hz for theta and 30-50Hz for gamma
%   or 'fast' 5-9Hz for theta and 70-100 Hz for gamma
%   also, input variable 'band' can be set up as a specificant passband
%   such as [59 61] as a notch filter for 60Hz noise removal. 
%
%
% David Wang 03/20/2019


function [amp_gamma,phase_theta]=PACprep(x,band,Fs)





if strcmp(band,'slow')
    gamma = [30 70]*2/Fs;
    theta = [2 5]*2/Fs;
else
    
    if strcmp(band,'fast')
        gamma = [70 100]*2/Fs;
        theta = [5 9]*2/Fs;
    else
        
        if  strcmp(band,'slowfast')
            gamma = [70 100]*2/Fs;
            theta = [2 5]*2/Fs;
        else
            
            if  strcmp(band,'fastslow')
                gamma = [30 70]*2/Fs;
                theta = [5 9]*2/Fs;
            else
                gamma = band(1,:)*2/Fs;
                theta = band(2,:)*2/Fs;
            end
        end
    end
end



or_theta =1;
or_gamma =1;
[b_gamma,a_gamma] = butter(or_gamma,gamma);
[b_theta,a_theta] = butter(or_theta,theta);

x_gamma = zeros(size(x));
x_theta = zeros(size(x));
z_gamma = zeros(size(x));
z_theta = zeros(size(x));
phase_theta = zeros(size(x));
amp_gamma = zeros(size(x));
for i = 1:size(x,1)
    x_gamma(i,:) = filter(b_gamma,a_gamma,x(i,:));
    x_theta(i,:) = filter(b_theta,a_theta,x(i,:));
end

for i = 1:size(x,1)
    z_gamma(i,:) = hilbert(x_gamma(i,:));
    z_theta(i,:) = hilbert(x_theta(i,:));
    %     phase_the = angle(z_theta(i,:));
    %     sec_phe = abs(hilbert(phase_the));
    %     phase_theta(i,:) = sec_phe;
    
    % phase_theta(i,:) = angle(z_theta(i,:));
    phase_theta(i,:) = angle(z_theta(i,:));   % phase of theta
    amp_gamma(i,:)  = abs(z_gamma(i,:));       % amplitdue of gamma
    
end



end
