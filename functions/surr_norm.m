% This function computes normoalizes xPAC using a surrogate procedure.
% [M_norm_t,PfPha,totalMI,M_raw] = surr_norm(amplitude_t,phase_t)
% amplitude_t is the input matrix of amplitudes and phase_t is the input 
% matrix of phases, where each row is an individual trial.
% M_norm_t is the vector of output normalized Modulation Index for trials. 
% PfPha is the vector of preferred phase of coupling for trials.
% totalMI is the cell matrix of surrogate data where each cell contains all
% (250) surrogate trials for each indiviual EEG trial. 
% M_raw is the vector of raw modulation index without
% normalization (surrogate).
% David Wang 03/20/2019



function [M_norm_t,PfPha,totalMI,M_raw] = surr_norm(amplitude_t,phase_t)

PfPha = zeros(1,size(amplitude_t,1));
M_norm_t = zeros(1,size(amplitude_t,1));
M_raw = zeros(1,size(amplitude_t,1));
totalMI = cell(1,size(amplitude_t,1));
for i = 1: size(amplitude_t,1)
    amplitude = amplitude_t(i,:);
    phase = phase_t(i,:);
    numpoints=length(amplitude);
    
    numsurrogate=250;
    minskip=50;
    maxskip=400;
    
    skip=ceil(numpoints.*rand(numsurrogate*2,1));
    skip(skip>maxskip)=[];
    skip(skip<minskip)=[];
    skip=skip(1:numsurrogate,1);
    surrogate_m=zeros(numsurrogate,1);
    
    
    z=amplitude.*exp(1i*phase);
    %% mean of z over time, prenormalized value
    m_raw=mean(z);
    
    %% compute surrogate values
    for s=1:numsurrogate
        surrogate_amplitude=[amplitude(skip(s):end) amplitude(1:skip(s)-1)];
        surrogate_m(s)=abs(mean(surrogate_amplitude.*exp(1i*phase)));
        %disp(numsurrogate-s)
    end
    %% fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
    % [surrogate_mean,surrogate_std]=normfit(surrogate_m);
    %% normalize length using surrogate data (z-score)
    % surrogate_mean = mean(surrogate_m);
    % surrogate_std = std(surrogate_m);
    [surrogate_mean,surrogate_std]=normfit(surrogate_m);
    m_norm_length=(abs(m_raw)-surrogate_mean)/surrogate_std;
    m_norm_phase=angle(m_raw);
    m_norm=m_norm_length*exp(1i*m_norm_phase);
    
    M_raw(i) = abs(m_raw);
    M_norm = abs(m_norm);
    PfPha(i) = angle(mean(z));
    M_norm_t(i) = M_norm;
    totalMI{i}=surrogate_m;
end
end
