%% By Mattia
function [stimulationTimeStamps] = SUAgetStimulationTimeStamps(experiment, save_data,Path)

%Path = get_path;
parameters=get_parameters;

filename = strcat(experiment.path, filesep, experiment.name, filesep, 'STIM1D.ncs');

FieldSelectionArray     = [0 0 1 0 1]; %TimeStamps, ChannelNumbers, SampleFrequencies, NumberValidSamples, Samples
ExtractHeaderValue      = 0;
ExtractMode             = 1;
ExtractModeArray        = [];
[SampleFrequencies, Samples] = Nlx2MatCSC(filename, FieldSelectionArray, ExtractHeaderValue, ExtractMode, ExtractModeArray);

clear TimeStamps
f = median(SampleFrequencies);
clear SampleFrequencies

% rearray and adjust
Samples = reshape(Samples, 1, size(Samples, 1) * size(Samples, 2));
Samples = Samples ./ 32.81; %adjust to microVolt; only for LFP
Samples = Samples > max(Samples) * parameters.preprocessing.nlx_load_Opto.digital2binary.Threshold;

%% find stimulus period
StimDStart = find(diff(Samples) == 1) + 1;
StimDEnd = find(diff(Samples) == - 1) - 1;
StimDshortinterval = find(diff(StimDStart) < f) + 1;
StimDStart(StimDshortinterval) = [];
StimDEnd(StimDshortinterval - 1) = [];

StimStart = StimDStart; % leave it in timestamps
StimEnd = StimDEnd; % leave it in timestamps

stimulationTimeStamps = [StimStart',StimEnd'];

%% Save
if save_data == 0
    return
else
    if ~exist(strcat(Path.output,filesep,'SUAstimulatioTimeStamps', filesep))
    mkdir(strcat(Path.output,filesep,'SUAstimulatioTimeStamps', filesep));
    end
    save(strcat(Path.output,filesep,'SUAstimulatioTimeStamps', filesep, experiment.name,'.mat'),'stimulationTimeStamps');
end
end