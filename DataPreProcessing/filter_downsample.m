% downsample file with previous low pass filter at nyquist frequency
% Sebastian Bitzenhofer 2016-09

function [TimeStamps, Samples, fs] = filter_downsample(TimeStamps, Samples, fs, DownsampleFactor)

%Input 
%TimeStamps: timepoints of samples as row vector in milliseconds
%Samples: samples as row vector in microvolt
%SamplingFrequency: samples recorded per second
%DownsampleFactor: factor to downsample by as number

%Output 
%TimeStamps: timepoints of samples as row vector in milliseconds
%Samples: samples as row vector in microvolt
%SamplingFrequency: samples recorded per second

if DownsampleFactor~=1
    high_cut=fs/DownsampleFactor/2;
    Samples=ZeroPhaseFilterZeroPadding(Samples,fs,[0 high_cut]);
    
    TimeStamps=TimeStamps(1:DownsampleFactor:end);
    Samples=Samples(1:DownsampleFactor:end);
    fs=fs/DownsampleFactor;
end
