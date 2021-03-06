function main_function_detectEvents_new(Experiment,Path,band,param,electrodes,flag_time_path)
%written by Sabine Gretenkord - last modified 25.10.2016
%event detection based on Nicole's functions
%'detection_discont_events' and 'plot_event_detection'

num_sd=param.th_r;
ExtractModeArray=param.ExtractModeArray;
fs=param.fs;
DownsampleFactor=param.DownsampleFactor;
% recalculate or not
FieldName=fieldnames (param);
if ismember('flag_recalculate',FieldName)
    flag_recalculate=param.flag_recalculate;
else
    flag_recalculate=1;
end

for iExperiment=1:length(Experiment)
    experiment=Experiment(iExperiment);
    iExperiment
    
    filename=Experiment(iExperiment).name;
    if isempty(electrodes)
        Electrodes=[experiment.HPreversal experiment.PL experiment.PL+12];
        Electrodes=[experiment.HPreversal experiment.PL];
    else
        Electrodes=[experiment.HPreversal electrodes];
    end
    
    FieldName=fieldnames (param);
    if ismember('flag_region',FieldName)
        if ismember('LEC',param.flag_region)
            csc5=experiment.LECProbereversal-8;
            csc2=experiment.LECProbereversal-11;
            Electrodes=[Electrodes csc5 csc2];
        end
    end
    
    numchannels=length(Electrodes);
    for ichan=1:numchannels
        chanChoice=Electrodes(ichan);
        if chanChoice~=0 % 0 means no channel in this region
            resultFolder=strcat(Path.temp,filesep,'detection_discont_events_fs_',num2str(fs),'_thr',num2str(num_sd),'\',filename,'\ch',num2str(chanChoice),'\');
            if exist (resultFolder) && ~flag_recalculate
                disp ('Exist, no recalculation')
            else
                mkdir(resultFolder)
                cd(resultFolder)
                
                %load signal, time vector and sampling rate
                File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(chanChoice),'.ncs');
                
                if nargin < 6
                    [T, Samples, fs] = load_nlx(File,ExtractModeArray);
                else
                    [T, Samples, fs] = load_nlx_stimulation_baseline(experiment, chanChoice,flag_time_path);
                end
                TimeStamps=T.TimeStamps;
                [TimeStamps, Samples, fs] = filter_downsample(TimeStamps, Samples, fs, DownsampleFactor);
                signal=ZeroPhaseFilter(Samples,fs,[band(1) band(2)]);
                time=TimeStamps;
                
                % cut signal or not
                if ismember('flag_cutSignal',FieldName)
                    signal=signal (1:param.flag_cutSignal*fs*60);
                    time =time(1:param.flag_cutSignal*fs*60);
                end
                
                if ~isempty(signal)
                    
                    %detect discontinuous events with Nicole's method
                    timeAnalyzed=length (time)/fs; % s
                    [osc_start,osc_end,oscLogi,rms,thr,minint,mindur,num_sd]  = detection_discont_events(signal,fs,num_sd); %samples
                    
                    %plot event detection
                    plot_event_detection(signal,time,fs,osc_start,osc_end,rms,thr);
                    title( strcat(  strrep(filename,'_','-') ,chanChoice  ) )
                    xlabel('Time(points)')
                    
                    %calculate parameters
                    oscDetect=[osc_start;osc_end];
                    silDetect=[];
                    for seg=1:length(osc_start)-1
                        silDetect(:,seg)=[osc_end(seg)+1;osc_start(seg+1)-1];
                    end
                    
                    eventFreq=size(oscDetect,2)/timeAnalyzed;
                    eventLength=mean(diff(oscDetect)/fs);
                    
                    %save detection result and figure
                    
                    save('oscDetect.mat','oscDetect','oscLogi','silDetect','timeAnalyzed','eventFreq','eventLength','minint','mindur','num_sd','thr')
                    saveas(gcf,'OscDetectFig')
                    close all
                    
                end
                
            end
            
            
        end
        %         clearvars -except flag_recalculate param Experiment DownsampleFactor iExperiment fs Path ichan Band num_sd electrodes ExtractModeArray band flag_time_path
    end
    
end











