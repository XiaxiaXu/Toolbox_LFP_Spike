function SUAgetStimulationProperties(experiments, save_data, repeatCalc,Path)
%Path = get_path;
for n_animal = 1:length(experiments)
    experiment = experiments(n_animal);
    if ~isempty(experiment.animal_ID)
        if  repeatCalc == 0 && exist(strcat(Path.output, filesep, 'StimulationPropertiesSUA', ...
                filesep, experiment.name, filesep, 'StimulationProperties_raw','.mat'))
            disp(['Already calculated stimProperties for: ' experiment.name ', expNumber ' num2str(n_animal)])
        else
            disp(['Calculating stimProperties for: ' experiment.name ', expNumber ' num2str(n_animal)])
            [StimTimestamps] = SUAgetStimulationTimeStamps(experiment, save_data,Path);
            StimTimestamps = StimTimestamps / 10;
            [~, signalD, ~] = nlx_load_Opto(experiment, 'Stim1D', [], 10, 0);
            signalD = digital2binary(signalD);
            [~, signalA, samplingrate] = nlx_load_Opto(experiment, 'Stim1A', [], 10, 0);
            
            stim_digital_single     = signalD;
            stim_digital_periode    = signalD;
            clearvars signalD
            
            %% Find single stimuli
            stimStart                   = find(diff(stim_digital_single) == 1);
            stimEnd                     = find(diff(stim_digital_single) == -1);
            stim_dur                    = (stimEnd-stimStart) / samplingrate;
            
            %% find stimulus period
            StimulationProperties_raw = {};
            
            %min interval (If interval is too small <1sec, make it longer
            stim_periode_Start        = find(diff(stim_digital_periode)==1);
            stim_periode_End          = find(diff(stim_digital_periode)==-1);
            stim_periode_interval     = (stim_periode_Start(2:end)-stim_periode_End  (1:end-1))/samplingrate;
            
            for k           = find(stim_periode_interval<1)
                stim_digital_periode(stim_periode_End(k) : stim_periode_Start(k+1)) = 1;
            end
            
            %min duration (If duration is too short <0.5sec, remove it)
            stim_periode_Start        = find(diff(stim_digital_periode)==1);
            stim_periode_End          = find(diff(stim_digital_periode)==-1);
            stim_periode_dur          = (stim_periode_End - stim_periode_Start)/samplingrate;
            
            for k           = find(stim_periode_dur<1)
                stim_digital_periode(stim_periode_Start(k):stim_periode_End  (k)) = 0;
            end
            
            stim_periode_Start        = find(diff(stim_digital_periode)==1);
            stim_periode_End          = find(diff(stim_digital_periode)==-1);
            
            
            %% Calculate period properties
            for j = 1 : length(stim_periode_Start)
                
                periodeStart    = stim_periode_Start(j);
                periodeEnd      = stim_periode_End(j);
                
                %% calculate n_pulses
                if j ==1
                    n_pulses   = max(find(stimStart < periodeEnd));
                else
                    n_pulses   = max(find(stimStart < periodeEnd)) - max(find(stimStart < periodeStart));
                end
                
                %% stim On times for each period
                
                if j ==1
                    singlePulseStartTime = (stimStart(1:max(find(stimStart < periodeEnd))));
                else
                    singlePulseStartTime = (stimStart(max(find(stimStart < periodeStart))+1:max(find(stimStart < periodeEnd))));
                end
                
                %% calculate mean pulse duration
                if j == 1
                    pulse_duration = round(mean(stim_dur(1:n_pulses))*1000); % give output in ms
                else
                    pulse_duration = round(mean(stim_dur(max(find(stimStart == periodeStart)):max(find(stimStart < periodeEnd))))*1000); % give output in ms
                end
                
                %% calculate frequency
                
                f_stim = round(n_pulses / round(stim_periode_dur(j)));
                
                %% calculate laser power
                laser_power = max(signalA(stim_periode_Start(j):stim_periode_End(j)));
                laser_power = round(-0.0001*(laser_power)^2 + 0.6228*laser_power - 3.2176); %% The formula is generated from measurements from laser
                
                %% identify stimulus type (frequency, constant, ramp, or sinus)
                % stim_class [0, 1, 2, 3, 4] = [Unknown, frequency, constant, ramp, sinus]
                if n_pulses > 1
                    % Is there more than 1 stimulation in the stimulation period?
                    stim_class = 'squareFreq';
                else
                    pulse_duration = NaN;
                    % calculate max/mean ratio, for constant light = 0.5, for ramp/sinus = 1
                    halfmax_analog = max(signalA(stim_periode_Start(j):stim_periode_End(j)))/2;
                    mean_analog    = mean(signalA(stim_periode_Start(j):stim_periode_End(j)));
                    max_mean_ratio = halfmax_analog/mean_analog;
                    if max_mean_ratio <= 0.55 % this number can be changed to fit data best. Ideally 0.50 but higher due to noise
                        stim_class = 'constant';
                        f_stim = NaN;
                    else
                        %% calculate number of zero cross going from - to + and + to -, if == then sinus, if not, something else.
                        signal_half = (signalA(stim_periode_Start(j):stim_periode_End(j))-signalA(stim_periode_Start(j))) - mean(signalA(stim_periode_Start(j):stim_periode_End(j))-signalA(stim_periode_Start(j)))*(5/3); % The *5/3 at end is added to prevent influence of noise
                        find_zero   = diff(sign(signal_half));
                        indx_up     = find(find_zero>0); %find all upward going zeros
                        indx_down   = find(find_zero<0); %find all downward going zeros
                        halfstimdur = round((stim_periode_End(j)-stim_periode_Start(j))/2);
                        maxfirsthalf = max(signalA(stim_periode_Start(j):stim_periode_End(j)-halfstimdur));
                        maxlasthalf = max(signalA(stim_periode_End(j)-halfstimdur:stim_periode_End(j)));
                        if round(maxfirsthalf) == round(maxlasthalf) && numel(indx_up) > 0
                            f_stim  = round(length(indx_up)/round(stim_periode_dur(j))); %correct frequency for sinus
                            indx_start = indx_up(2)-indx_up(1);
                            indx_end = indx_up(end-1)-indx_up(end-2);
                            if abs(indx_start-indx_end)<100  % check if difference between 2 first and 2 last indx is equal or not (equal = sinus, unequal = chirp)
                                stim_class = 'sinus';
                            else
                                stim_class = 'chirp';
                            end
                        else
                            % Check if first number is smaller than last
                            if signalA(stim_periode_Start(j)) < signalA(stim_periode_End(j))
                                stim_class = 'ramp';
                                f_stim = NaN;
                            else
                                stim_class = 0;
                            end
                        end
                    end
                end
                
                %% store in file
                stimPeriodeStart = StimTimestamps(j,1);
                stimPeriodeEnd = StimTimestamps(j,2);
                StimulationProperties_raw(j,:) = ...
                    {stimPeriodeStart, stimPeriodeEnd,round(stim_periode_dur(j)), n_pulses, pulse_duration, f_stim, laser_power, stim_class,singlePulseStartTime};
            end
            %% Save
            if ~isempty(StimulationProperties_raw)
                if save_data == 0
                elseif save_data == 1
                    if ~exist(strcat(Path.output,filesep,'StimulationPropertiesSUA', filesep, experiment.name, filesep))
                        mkdir(strcat(Path.output,filesep,filesep,'StimulationPropertiesSUA', filesep, experiment.name, filesep));
                    end
                    save(strcat(Path.output,  filesep, 'StimulationPropertiesSUA', filesep, experiment.name, filesep, 'StimulationProperties_raw','.mat'),'StimulationProperties_raw');
                end
            else
                disp('StimProperties was not calculated')
            end
        end
    end
    
end
end