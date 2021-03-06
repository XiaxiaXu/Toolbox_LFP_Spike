function [T, signal, samplingrate] =load_nlx_stimulation_baseline(experiment, CSC, flag_time_path)
% full signal
downsample_factor=1;
load(strcat(flag_time_path, filesep, experiment.name, filesep, 'StimulationProperties_corrected'));

[time, ~, ~] = nlx_load_Joachim(experiment,'stim1D',[],10,0);

T_first = time(1);
T_last = time(end);
clearvars time

P_start = [T_first cell2mat(StimulationProperties_corrected(:,2))' T_last];
P_baseline = find(diff(P_start) > (15*60*10^6));
P_baselines = [P_start(P_baseline(1)) P_start(P_baseline(end))];

 bb = 1;
        ExtractModeArray(1) = P_baselines(bb);
        ExtractModeArray(2) = P_baselines(bb)+15*60*10^6;
        [timeX, signalX,samplingrate] = nlx_load_Joachim(experiment,CSC,ExtractModeArray,downsample_factor,0);
        timeX = timeX(1:(15*60*samplingrate));
        signalX = signalX(1:(15*60*samplingrate));
        signal = signalX(samplingrate*10:end-samplingrate*10); % cut away 10 sec from start and end. This is to prevent artifact/filter artifacts coming from power estimation or end of last trial.
        time = timeX(samplingrate*10:end-samplingrate*10);
        
    fs=samplingrate;
    TimeStamps= 10^3/fs:10^3/fs:(length(signal)*10^3/fs);% ms, no matter when the recording started, just make the first start point as 1/fs
    TimeStamps_true=time;
    
    T.TimeStamps=TimeStamps; 
    T.TimeStamps_true=TimeStamps_true; % s
