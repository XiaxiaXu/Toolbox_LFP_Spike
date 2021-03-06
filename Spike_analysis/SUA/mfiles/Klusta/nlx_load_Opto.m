%% By Joachim Ahlbeck & Nicola improved by Sebastian Bitzehnofer
function [time, signal, fs, high_cut] = nlx_load_Opto(experiment, CSC, ExtractModeArray, downsampling_factor, save_data)
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

parameters  = get_parameters;
Path        = get_path;

%% Script
if isa(CSC,'numeric')
     filename = strcat(experiment.path,experiment.name, '\CSC',num2str(CSC),'.ncs');
%      filename = strcat(experiment.path, '\CSC',num2str(CSC),'.ncs');
elseif isa(CSC,'char')
    filename = strcat(experiment.path,experiment.name,filesep,CSC,'.ncs');
end
    
FieldSelectionArray     = [1 0 0 0 1]; %TimeStamps, ChannelNumbers, SampleFrequencies, NumberValidSamples, Samples
ExtractHeaderValue      = 1;
if isempty(ExtractModeArray)
    ExtractMode         = 1;
else
    ExtractMode         = 2;
end
[time, signal, Header] = Nlx2MatCSC(filename, FieldSelectionArray, ExtractHeaderValue, ExtractMode, ExtractModeArray);

%rearray and adjust
signal                = reshape(signal,1,size(signal,1)*size(signal,2));
signal                = signal./parameters.preprocessing.nlx_load_Opto.voltAdjust; % Adjust to microVolt, normally 32.81

%time and samplingrate
fs=str2double(Header{14,1}(20:end));
if isnan(fs)
    fs=str2double(Header{15,1}(20:end));
end    
time      = linspace(time(1),time(end)+512/fs*10^6,length(signal)); %s

% low pass filter at nyquist frequency and downsample
high_cut=[];
if downsampling_factor~=1
    if isa(CSC,'numeric')
        high_cut = fs/downsampling_factor/2;
        signal   = ZeroPhaseFilter(signal,fs,[0 high_cut]);
    end
    time     = time(1:downsampling_factor:end);                           %resample with resamp factor as an input instead of every 10th value INTRODUCES ERRORS
    signal   = signal(1:downsampling_factor:end);
    fs       = fs/downsampling_factor;
end

%% save
if save_data==1
    if ~exist(strcat(Path.temp,filesep,'nlx_load_Opto',filesep,experiment.name))
        mkdir(strcat(Path.temp,filesep,'nlx_load_Opto',filesep,experiment.name))
    end
    if isa(CSC,'numeric')
        savename=strcat(Path.temp,filesep,'nlx_load_Opto',filesep,experiment.name,filesep,'CSC',num2str(CSC));
    elseif isa(CSC,'char')
        savename=strcat(Path.temp,filesep,'nlx_load_Opto',filesep,experiment.name,filesep,CSC);
    end
    save(savename,...
        'time','signal','fs','high_cut', '-v7.3')
end
end