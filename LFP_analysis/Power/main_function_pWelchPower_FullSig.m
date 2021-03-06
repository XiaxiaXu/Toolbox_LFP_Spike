function main_function_pWelchPower_FullSig(Experiment,Path,params,electrodes,flag_time_path)
%[Power_FullSig,f,~] = pWelchSpectrum(LFP,params.win_sec,0,2024,fs,100);

ExtractModeArray=params.ExtractModeArray;
band=params.filter_band;
flag_CutData=params.flag_CutData;
DownsampleFactor=params.DownsampleFactor;
FieldName=fieldnames (params);

for iExperiment=1:length(Experiment)
    experiment=Experiment(iExperiment);
    iExperiment
    filename=Experiment(iExperiment).name;
    if isempty(electrodes)
        Electrodes=[experiment.HPreversal experiment.PL];
    else
        Electrodes=electrodes;
    end
    
     if ismember('flag_LEC',FieldName)
            csc=experiment.LECreversal-5;
            Electrodes=[Electrodes csc];
    end
    
    
    fs=params.fs;
    resultFolder=strcat(Path.output,filesep,'pWelchPower\',filename,'_fs',num2str(fs),filesep);
    
    numchannels=length(Electrodes);
    for ichan=1:numchannels
        chanChoice=Electrodes(ichan);
        CSC=chanChoice;
        %load signal, time vector and sampling rate
        File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(chanChoice),'.ncs');
        if nargin < 5
            [T, Samples, Fs] = load_nlx(File,ExtractModeArray);
        else
            [T, Samples, Fs] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
        end
        TimeStamps=T.TimeStamps;
        [~, samples, fs] = filter_downsample(TimeStamps, Samples, Fs, DownsampleFactor);
        signal=ZeroPhaseFilter(samples,fs,[band(1) band(2)]);
        
        if flag_CutData
            if length(signal)>flag_CutData(2)*60*fs
                signal=signal(flag_CutData(1)*60*fs:flag_CutData(2)*60*fs);
            else
                signal=signal(flag_CutData(1)*60*fs:end);
            end
        end
        
        [Power,f,~] = pWelchSpectrum(signal,params.win_sec,0,2024,fs,100);
         win_sec=  params.win_sec;    
        save_name=strcat('CSC',num2str(CSC));
         
        if ~exist( resultFolder )
            mkdir(resultFolder)
        end
        cd(resultFolder)
        save(strcat(save_name,'.mat'),'Power','f','win_sec','fs')
    end
end



        

