function SpikesMatrix = getSUAMatrix(Experiment, Path,StimRegion, RespArea)
%% By Xiaxia

SpikesMatrix=[];
for n_animal =1: length(Experiment)
    n_animal
    experiment = Experiment(n_animal);
    if nnz( ~ isnan(experiment.StimRegion))
        if ~isempty(experiment.animal_ID)
            if strcmp(RespArea, 'PFC2')
                foldersave = strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAPFC23',filesep,experiment.name);
            elseif strcmp(RespArea, 'PFC5')
                foldersave = strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAPFC56',filesep,experiment.name);
            elseif strcmp(RespArea, 'HP')
                foldersave = strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAHP',filesep,experiment.name);
            end
            
            if exist(strcat(foldersave,'\SUAdata_pulses.mat'))
                load(strcat(foldersave,'\SUAdata_pulses.mat'))
                spikes_matrix = SUAdata_pulses.pulse_spike_matrix;
                if size(spikes_matrix, 2) > 1
                    spikes_units = squeeze(mean(spikes_matrix));
                else
                    spikes_units = squeeze(mean(spikes_matrix))';
                end
                SpikesMatrix=[SpikesMatrix;spikes_units];
            end
        end
    end
end

%spy(SpikesMatrix,'.',6)



