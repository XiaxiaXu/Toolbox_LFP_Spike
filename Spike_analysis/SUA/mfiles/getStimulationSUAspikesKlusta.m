function SUAdata = getStimulationSUAspikesKlusta(experiment, save_data, repeatCalc,Path)
% by Mattia

% Path = get_path;

if exist(strcat(Path.output, filesep, 'StimulationSUAspikeTimesKlusta', filesep, ...
        experiment.name, filesep, 'StimulationSUAspikeTimes.mat')) && repeatCalc == 0
    load(strcat(Path.output, filesep, 'StimulationSUAspikeTimesKlusta', filesep, ...
        experiment.name, filesep, 'StimulationSUAspikeTimes.mat'))
else
    display(strcat('calculating spiketimes for: ', experiment.name))
    %     Path = get_path;
    if exist(strcat(Path.SUA,'\SortedSUA\', experiment.animal_ID, '.mat'))
        load(strcat(Path.SUA,'\SortedSUA\', experiment.animal_ID, '.mat'))
        
        if ~ isempty(fieldnames(SUAinfo{1}))
%                     for recording = 1 : numel(SUAinfo)
%                         if strcmp(SUAinfo{recording}(1).file, experiment.name)
%                             break
%                         end
%                     end
            recording = 1;
            SUA = SUAinfo{recording};
            load(strcat(Path.output, filesep, 'StimulationPropertiesSUA', filesep, experiment.name,  filesep, 'StimulationProperties_raw.mat'));
            StimulationProperties_corrected=StimulationProperties_raw;
            if ~ isempty(StimulationProperties_corrected)
                StimProp = round(cell2mat(StimulationProperties_corrected(:, 1 : 2))); % to ms
                first_stim = StimProp(1, 1); % dirty solution, doesn't generalize over non VIR mice, needs to be improved
                for idx = 1 :  size(SUA, 2)
                    if ~ isnan(SUA(idx).Timestamps)
                        try
                            spikes(idx, round(SUA(idx).Timestamps / 32)) = 1; % swithc to ms
                        catch
                            SUA(idx).Timestamps(SUA(idx).Timestamps == 0) = 1; % to avoid problems with indexing
                            spikes(idx, ceil(SUA(idx).Timestamps / 32)) = 1; % swithc to ms
                        end
                        channels(idx) = SUA(idx).channel;
                        clusters(idx) = SUA(idx).ClusterID;
                        if  first_stim >= 15*60*1000 %15min baseline
                            spikes_baseline = nnz(SUA(idx).Timestamps / 32 < 15*60*1000);
                            FiringRatebaseline(idx) = spikes_baseline / (15*60*1000 / 1000); % in Hz
                        else
                            FiringRatebaseline(idx) = NaN;
                        end
                    else
                        spikes(idx, :) = 0;
                    end
                end
            else
                StimulationProperties_corrected = [];
                idx = 0;
            end
            
            SUAdata.FiringRatebaseline=FiringRatebaseline;
            
            %            calculate the stimulation
            for i=1: length(StimulationProperties_corrected)
                
                if isnumeric(StimulationProperties_corrected{i,8})~=1
                    if  ~isempty(StimulationProperties_corrected{i,8})
                        stimTypes(i)=cellstr(StimulationProperties_corrected(i,8));
                    else
                        stimTypes(i)={'NaN'};
                    end
                end
            end
            Logi_ramp=find(strcmp(stimTypes,'ramp')==1);
            ramps=Logi_ramp; % count number of ramps
            ramp_spike_matrix = zeros(size(ramps,2), size(SUA, 2), 9000); % n_ramps x n_units x time %9s
            %         experiment.ramp = find(Logi_ramp==1);
            idx=0;
            if numel(ramp_spike_matrix) > 0 && size(ramps,2) > 1
                for stimulation =ramps
                    idx = idx + 1;
                    stim_end_ramp = round(cell2mat(StimulationProperties_corrected(stimulation, 2))); % stimproperties are with fs = 32000, hence the correction
                    if stim_end_ramp + 3000 > length(spikes)
                        stim_end_ramp = length(spikes) - 3000;
                    end
                    if ~ isempty(stim_end_ramp)
                        ramp_spike_matrix(idx, :, :) = spikes(:, stim_end_ramp-5999 : stim_end_ramp + 3000);
                    else
                        idx = idx - 1;
                    end
                end
                ramppre = sum(sum(ramp_spike_matrix(:, :, 1 : 2995), 3)); % light artifact
                rampduring = sum(sum(ramp_spike_matrix(:, :, 3005 : 5995), 3)); % light artifact
                ramppre_single_ramp = sum(ramp_spike_matrix(:, :, 1 : 2995), 3); % light artifact
                rampduring_single_ramp = sum(ramp_spike_matrix(:, :, 3005 : 5995), 3); % light artifact
                ramppvalue = zeros(1, length(ramppre));
                for unit = 1 : length(ramppre)
                    ramppvalue(unit) = signrank(ramppre_single_ramp(:, unit), ...
                        rampduring_single_ramp(:, unit));
                end
                OMI = (rampduring - ramppre)./ (rampduring + ramppre);
                OMI2= (rampduring - ramppre)./ ramppre;
            end
            Logi_square1=find(strcmp(stimTypes,'squareFreq')==1);
            pulses=Logi_square1;% count number of pulses
            pulse_spike_matrix = zeros(size(pulses,2),size(SUA, 2), 9000); % n_ramps x n_units x time %9s
            %         experiment.pulse = find(Logi_pulse==1);
            idx=0;
            if numel(pulse_spike_matrix) > 0 && size(pulses,2) > 1
                for stimulation = pulses
                    idx = idx + 1;
                    stim_end_pulse = round(cell2mat(StimulationProperties_corrected(stimulation, 2))); % stimproperties are with fs = 32000, hence the correction
                    if stim_end_pulse + 3000 > length(spikes)
                        stim_end_pulse = length(spikes) - 3000;
                    end
                    if ~ isempty(stim_end_pulse)
                        pulse_spike_matrix(idx, :, :) = spikes(:, stim_end_pulse - 5999 : stim_end_pulse + 3000);
                    else
                        idx = idx - 1;
                    end
                end
                pulsepre = sum(sum(pulse_spike_matrix(:, :, 1 : 2995), 3)); % light artifact
                pulseduring = sum(sum(pulse_spike_matrix(:, :, 3005 : 5995), 3)); % light artifact
                pulsepre_single_pulse = sum(pulse_spike_matrix(:, :, 1 : 2995), 3); % light artifact
                pulseduring_single_pulse = sum(pulse_spike_matrix(:, :, 3005 : 5995), 3); % light artifact
                pulsepvalue = zeros(1, length(pulsepre));
                for unit = 1 : length(pulsepre)
                    pulsepvalue(unit) = signrank(pulsepre_single_pulse(:, unit), ...
                        pulseduring_single_pulse(:, unit));
                end
                pulseOMI = (pulseduring - pulsepre)./ (pulseduring + pulsepre);
                pulseOMI2= (pulseduring - pulsepre)./ pulsepre;
                %
                SUAdata.ramp_spike_matrix = ramp_spike_matrix;
                SUAdata.channels = channels;
                SUAdata.clusters = clusters;
                SUAdata.rampOMI = OMI;
                SUAdata.rampOMI2 = OMI2;
                SUAdata.ramppvalue=ramppvalue;
                SUAdata.pulse_spike_matrix = pulse_spike_matrix;
                SUAdata.pulseOMI = pulseOMI;
                SUAdata.pulseOMI2 = pulseOMI2;
                SUAdata.pulsepvalue = pulsepvalue;
            else
                SUAdata.ramp_spike_matrix = [];
                SUAdata.pulse_spike_matrix = [];
                SUAdata.channels = [];
                SUAdata.clusters = [];
                SUAdata.rampOMI = [];
                SUAdata.rampOMI2= [];
                SUAdata.ramppvalue=[];
                SUAdata.pulseOMI = [];
                SUAdata.pulseOMI2= [];
                SUAdata.FiringRatebaseline = [];
                SUAdata.pulsepvalue = [];
            end
        end
                
        if save_data == 0
            disp('DATA NOT SAVED!');
        elseif save_data==1
            if ~exist(strcat(Path.output,filesep,'StimulationSUAspikeTimesKlusta',filesep,experiment.name))
                mkdir(strcat(Path.output,filesep,'StimulationSUAspikeTimesKlusta',filesep,experiment.name));
            end
            save(strcat(Path.output, filesep, 'StimulationSUAspikeTimesKlusta', filesep, ...
                experiment.name, filesep, 'StimulationSUAspikeTimes'), 'SUAdata')
        end
        
    else
        SUAdata=[];
    end
        
end

