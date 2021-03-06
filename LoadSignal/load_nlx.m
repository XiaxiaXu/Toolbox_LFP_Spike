% Load neuralynx file (.nlx) to matlab (.mat)
% Sebastian Bitzenhofer 2016-09

function [T, Samples, fs] = load_nlx(File,ExtractModeArray)

%Input 
%file: file path and name as string
%ExtractModeArray: []-full trace; [start_time end_time] in milliseconds

%Output 
%TimeStamps: TimeStamps of samples as row vector in
%milliseconds [1xlength samples]
%Samples: samples as row vector in microvolt [1xlength samples]
%fs: samples recorded per second [1x1] 

FieldSelectionArray     = [1 0 0 0 1]; %TimeStamps, ChannelNumbers, SampleFrequencies, NumberValidSamples, Samples
ExtractHeaderValue      = 1;
if isempty(ExtractModeArray)
    ExtractMode         = 1;
    [TimeStamps, Samples, Header] = Nlx2MatCSC(File, FieldSelectionArray, ExtractHeaderValue, ExtractMode, ExtractModeArray);
else
    ExtractMode         = 4;
    ExtractModeArrayCorr(1)=(ExtractModeArray(1))*10^3;
    ExtractModeArrayCorr(2)=(ExtractModeArray(2))*10^3;
    [TimeStamps, Samples, Header] = Nlx2MatCSC(File, FieldSelectionArray, ExtractHeaderValue, ExtractMode, ExtractModeArrayCorr);
end

%rearray and adjust
ADBitVolts=str2double(Header{~cellfun('isempty', strfind(Header,'-ADBitVolts '))}(length('-ADBitVolts '):end));
fs=str2double(Header{~cellfun('isempty', strfind(Header,'-SamplingFrequency '))}(length('-SamplingFrequency '):end));

Samples=reshape(Samples,1,size(Samples,1)*size(Samples,2)).*ADBitVolts*10^6; % Adjust to microVolt
TimeStamps_true=linspace(TimeStamps(1),TimeStamps(end)+511/fs*10^6,length(Samples))/10^3; %adjust to msec; very small errors indtroduced due to imperfect sampling of the neuralynx system (examples tested with ~1 wrong timestamp in 1.5 hour recording leading to 2 ms difference in total duration)
%TimeStamps_true=linspace(TimeStamps(1),TimeStamps(end)+511/fs,length(Samples));%s; very small errors indtroduced due to imperfect sampling of the neuralynx system (examples tested with ~1 wrong timestamp in 1.5 hour recording leading to 2 ms difference in total duration)
TimeStamps=10^3/fs:10^3/fs:length(Samples)*10^3/fs;% ms, no matter when the recording started, just make the first start point as 1/fs

%correct for unprecise loading of Nlx2Mat with ExtractModeArray
if ~isempty(ExtractModeArray) % for example (15*60*10^6)
    timecorrect1=find(TimeStamps_true>=ExtractModeArray(1),1);
    timecorrect2=find(TimeStamps_true>=ExtractModeArray(2),1);
    
    Samples=Samples(timecorrect1:timecorrect2);
    TimeStamps_true=TimeStamps_true(timecorrect1:timecorrect2);
    TimeStamps=10^3/fs:10^3/fs:length(Samples)*10^3/fs; % ms
end

T.TimeStamps_true=TimeStamps_true*10^3; % divide 10^6 to second
T.TimeStamps=TimeStamps; % ms
%T.TimeStamps_true_Start=TimeStamps_true(1)*10^3; % mill second, divide 10^6 to second
