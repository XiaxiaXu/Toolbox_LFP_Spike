%% By Joachim Ahlbeck & Nicola improved by Sebastian Bitzehnofer
function [time, signal, fs] = nlx_load_real_time(experiment, CSC, ExtractModeArray, downsampling_factor)
%%%%%%%%%%%%%%%% Convertion of recording/stimulus channel %%%%%%%%%%%% 
%import nlx file and convert to matlab format
% This function require: Nlx2MatCSC from neuralynx and ZeroPhaseFilter
    % Can be found here: Q:\RG Opatz Super Analyser Pro 2015 (ROSA)\ExternalToolboxes\Neuralynx to Matlab Export
% Example of inputs:
    % number of experiment in get_experiment_list
    % CSC: Numeric channel numbers for CSC/ channel name as string
    % ExtractModeArray: []-full trace; [start_time end_time]-only loads time window (not precise since only nlx timestamps used: load little broader window and readjust)
    % downsampling_factor (e.g. 10) (1 to not downsample)
    % save_data: 1-save, 0-don't
% Example: [time, signal, fs] = nlx_load_OptoPFC(1, 1, [], 10, 0)
% Example: [time, signal, fs] = nlx_load_OptoPFC(1, 'STIM1D', [], 1, 0)

%%
if isa(CSC,'numeric')
     filename = strcat(experiment.path,filesep, experiment.name, '\CSC',num2str(CSC),'.ncs');
%      filename = strcat(experiment.path, '\CSC',num2str(CSC),'.ncs');
elseif isa(CSC,'char')
    filename = strcat(experiment.path,filesep, experiment.name, filesep, CSC,'.ncs');
end

FieldSelectionArray     = [1 0 0 0 1]; %TimeStamps, ChannelNumbers, SampleFrequencies, NumberValidSamples, Samples
ExtractHeaderValue      = 1;
if isempty(ExtractModeArray)
    ExtractMode         = 1;
    [time, signal, Header] = Nlx2MatCSC(filename, FieldSelectionArray, ExtractHeaderValue, ExtractMode, ExtractModeArray);
else
    ExtractMode         = 4;
    ExtractModeArrayCorr(1)=(ExtractModeArray(1)-100)*10^3;
    ExtractModeArrayCorr(2)=(ExtractModeArray(2)+100)*10^3;
    [time, signal, Header] = Nlx2MatCSC(filename, FieldSelectionArray, ExtractHeaderValue, ExtractMode, ExtractModeArrayCorr);    
end

%rearray and adjust
ADBitVolts=str2double(Header{~cellfun('isempty', strfind(Header,'-ADBitVolts '))}(length('-ADBitVolts '):end));
fs=str2double(Header{~cellfun('isempty', strfind(Header,'-SamplingFrequency '))}(length('-SamplingFrequency '):end));

signal=reshape(signal,1,size(signal,1)*size(signal,2)).*ADBitVolts*10^6; % Adjust to microVolt
time=linspace(time(1),time(end)+511/fs*10^6,length(signal)); %adjust to s; 

%correct for unprecise loading of Nlx2Mat with ExtractModeArray
if ~isempty(ExtractModeArray)
    timecorrect1=find(time>=ExtractModeArray(1),1);
    timecorrect2=find(time>=ExtractModeArray(2),1);
    
    time=time(timecorrect1:timecorrect2);
    signal=signal(timecorrect1:timecorrect2);
end

% low pass filter at nyquist frequency and downsample
if downsampling_factor~=1
    if isa(CSC,'numeric')
        high_cut = fs/downsampling_factor/2;
        signal   = ZeroPhaseFilter(signal,fs,[0 high_cut]);
    end
    time     = time(1:downsampling_factor:end);                           %resample with resamp factor as an input instead of every 10th value INTRODUCES ERRORS
    signal   = signal(1:downsampling_factor:end);
    fs       = fs/downsampling_factor;
end
