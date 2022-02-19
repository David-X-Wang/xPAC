function [MI_t,PfPha]=ModIdxRT(phase_t,amplitude_t,nBins)
%
% Function to compute Modulation Index (measure of Phase Amplitude Coupling)
% - Ideally extract the phase and amplitude using hilbert
% - Based on Tort et.al 2010, uses Kullback-Leibner distance to measure
% disparity between two distribution
% - nBin defaults to 18 bins spaced between -180 to 180 (-pi to pi).
% changing bins
%
% Ryan Tan 10/22/18
if nargin<3
    nBins=18; % default value
end

MI_t = zeros(1,size(amplitude_t,1));
distKL_t = zeros(1,size(amplitude_t,1));
PfPha = zeros(1,size(amplitude_t,1));
for i = 1: size(amplitude_t,1)
    ampl = amplitude_t(i,:);
    phas = phase_t(i,:);
    
    %% bin the amplitudes according to the phases
    
    binEdges=linspace(-pi,pi,nBins+1);
    %binCenters=binEdges(1:end-1)-diff(binEdges)/2;
    
    [~,binIdx]=histc(phas,binEdges);
    
    ampBin=zeros(1,nBins);
    for bin=1:nBins
        if any(binIdx==bin)
            ampBin(bin)=mean(ampl(binIdx==bin));
        end
    end
    
    ampSig=ampBin/sum(ampBin);
    amplNorm=ones(1,nBins)./nBins;
    
    % prefered phase bin (where has strongest MI)

    Phaspan = -pi:2*pi/(nBins-1):pi;
    phasecandidate = Phaspan(ampBin==max(ampBin));
    if length(phasecandidate)~=1
        continue
    else
        PfPha(i) = phasecandidate;
    end
    % in the special case where observed probability in a bin is 0, this tweak
    % allows computing a meaningful KL distance nonetheless
    if any(ampSig==0)
        ampSig(ampSig==0)=eps;
    end
    
    %% KL Distance and ML
    distKL=sum(ampSig.*log(ampSig./amplNorm));
    
    MI=distKL./log(nBins);
    MI_t(i) = MI;
    distKL_t(i) = distKL;
end
end

