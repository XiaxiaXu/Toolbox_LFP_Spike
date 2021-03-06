function main_function_get_oscillation_properties(Experiment,Path,params,electrodes,flag_time_path)

DownsampleFactor=params.DownsampleFactor;
num_sd=params.th_LFP; % 2~2.5
band=params.band;
flag_band=params.flag_band;  % '_gamma';
ExtractModeArray=params.ExtractModeArray;
flag_CutData=params.flag_CutData;
FieldName=fieldnames (params);

for iExperiment=1:length(Experiment)
    experiment=Experiment(iExperiment);
    iExperiment
    filename=Experiment(iExperiment).name;
    
    if isempty(electrodes)
        Electrodes=[experiment.HPreversal experiment.PL];
    else
        Electrodes=[experiment.HPreversal electrodes];
    end
    
    if ismember('flag_LEC',FieldName)
            csc=experiment.LECProbereversal-12:experiment.LECProbereversal-5;
            %Electrodes=[Electrodes csc];
            Electrodes=[csc];
    end
    
    numchannels=length(Electrodes);
    
    for ichan=1:numchannels
        CSC=Electrodes(ichan);
        %load signal, time vector and sampling rate
        File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(CSC),'.ncs');
        
        if nargin < 5
            [T, recordingRaw,samplingrate_MUA] = load_nlx(File,ExtractModeArray);
            TimeStamps=T.TimeStamps;
            [time, Samples, fs] = filter_downsample(TimeStamps, recordingRaw, samplingrate_MUA, DownsampleFactor);
            signal=ZeroPhaseFilter(Samples,fs,band);
            
            if flag_CutData
                Start_Point=60*fs;
                End_Point=14*60*fs;
                try
                    time=time(Start_Point:End_Point);
                    signal=signal(Start_Point:End_Point);
                catch
                    time=time(Start_Point:end);
                    signal=signal(Start_Point:end);
                end
            end
            
        else
            [T, recordingRaw,samplingrate_MUA] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
            TimeStamps=T.TimeStamps;
            [time, Samples, fs] = filter_downsample(TimeStamps, recordingRaw, samplingrate_MUA, DownsampleFactor);
            signal=ZeroPhaseFilter(Samples,fs,band);
        end
        
        timeAnalyzed=time(end)-time(1);
        [oscStart,oscEnd,~,~,~,~,~,~]  = detection_discont_events(signal,fs,num_sd); %samples
        
        %% get OscProperties
        [OscStructure] = getOscInformation_baseline(time,signal, fs,oscStart,oscEnd); % get occurence, amplitude and duration of oscillations
        OscOccurrence=OscStructure.occurrenceOsc; % store vars for each animal for laterplotting
        OscAmplitude=OscStructure.meanAmplOsc(1);
        OscDuration=OscStructure.meanDurOsc(1);
        
        % save the results
        if ~exist(  strcat(strcat(Path.output,filesep,'OscProperties',filesep,filename))  )
            mkdir(strcat(Path.output,filesep,'OscProperties',filesep,filename))
        end
        % mkdir(strcat(Path.output,filesep,'Firing_rate',filesep,filename))
        cd(strcat(Path.output,filesep,'OscProperties',filesep,filename))
        save(strcat('CSC',num2str(CSC),flag_band),'OscOccurrence','OscAmplitude','OscDuration')
    end
end


