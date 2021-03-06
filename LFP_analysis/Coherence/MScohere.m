function [Cxy,freq] = MScohere(signal1,signal2,windowSize,overlap,nfft,fs)            
%% Calculate Magnitude squared coherence using matlab function mscohere
%       By Joachim
%       This function can be used both with continous and discontinous
%       activity. Note: for disc. act. overlap must be = 0.

%       INPUTS:
%       signal1,signal2: signal in row vector
%       windowSize:      windowSize in sec (2)
%       overlap:         overlap in sec    (0)
%       nfft:            FFT lenght        (2048
%       fs:              samplingfrequency
[Cxy, freq] = mscohere(signal1, signal2, hanning(windowSize*fs), overlap*fs, nfft, fs);
end