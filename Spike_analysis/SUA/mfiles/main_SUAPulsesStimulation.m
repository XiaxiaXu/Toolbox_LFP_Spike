function main_SUAPulsesStimulation(Experiment,Path,StimRegion,RespArea,save_data, repeatCalc)

% compute all sort of computations, on an experiment by experiment basis
for exp_idx = 1: size(Experiment, 2)
    exp_idx
    % select experiment
    experiment = Experiment(exp_idx);
    if strcmp(RespArea, 'PFC5') % PFC2 PFC5 HP PFC;
        folder2save = strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAPFC56',filesep,experiment.name);
    elseif strcmp(RespArea, 'PFC2')
        folder2save =strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAPFC23',filesep,experiment.name);
    elseif strcmp(RespArea, 'HP')
        folder2save =strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAHP',filesep,experiment.name);
    end
    folder4stim_props =strcat(Path.output, filesep, 'StimulationPropertiesSUA', filesep, experiment.name);
   
    % compute power for ramps
    %     disp(['Computing Ramp Power stuff for animal ', experiment.animal_ID])
    %     for CSC = str2num(experiment.CSC)
    %         getStimulationPowerSingleRamps(experiment, CSC, save_data, repeatCalc,folder4stim_props, folder4power_ramps)
    %     end
    % extract spike matrices for ramps
    disp(['Extracting spike matrices for animal ', experiment.animal_ID])
    %     SUAdata_ramps = getRampSpikeMatrix(experiment, save_data, repeatCalc, folder4ramps, folderSM, folder4stim_props);
    SUAdata_pulses = getPulseSpikeMatrix(experiment, save_data, repeatCalc,  folder2save, folder4stim_props,RespArea,Path);
end
