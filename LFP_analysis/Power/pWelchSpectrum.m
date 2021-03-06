function [Pxx,freq,Pxxc] = pWelchSpectrum(signal,windowSize,overlap,nfft,fs,maxFreq)
%% Calculate powerSpectrum using matlab function pWelch
%       By Joachim
%       This function can be used both with continous and discontinous
%       activity. Note: for disc. act. overlap must be = 0.

%       INPUTS:
%       signal: ´        signal in row vector
%       windowSize:      windowSize in sec (2)
%       overlap:         overlap in sec    (0)
%       nfft:            FFT lenght        (2048)
%       fs:              samplingfrequency
[Pxx,freq,Pxxc] = pwelch(signal, hanning(windowSize*fs), overlap*fs, nfft, fs,'ConfidenceLevel',0.95);
freq=freq(freq<=maxFreq)';
Pxx=Pxx(1:length(freq))';
Pxxc=Pxxc(1:length(freq),:)';
end