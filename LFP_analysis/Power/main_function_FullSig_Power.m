function main_function_FullSig_Power(Experiment,Path,params,electrodes,flag_time_path)

fs=params.fs;
ExtractModeArray=params.ExtractModeArray;
DownsampleFactor=params.DownsampleFactor;
flag_CutData=params.flag_CutData;

FieldName=fieldnames (params);
if ismember('flag_recalculate',FieldName)
    flag_recalculate=params.flag_recalculate;
else
    flag_recalculate=1;
end

for iExperiment=1:length(Experiment)
    iExperiment
    filename=Experiment(iExperiment).name;
    experiment=Experiment(iExperiment);
    
    FieldName=fieldnames (params);
    
    % set channel
    if isempty(electrodes)
        
        if ismember('flag_region',FieldName)
            if ismember('LEC5',params.flag_region)
                csc=experiment.LECProbereversal-8;
            end
            if ismember('LEC2',params.flag_region)
                csc=experiment.LECProbereversal-11;
            end
            if ismember('vCA1',params.flag_region)
                csc=experiment.HPreversal;
            end
            if ismember('PFC2',params.flag_region)
                csc=experiment.PL;
            end
            if ismember('PFC5',params.flag_region)
                csc=experiment.PL+12;
            end
        else
            csc=[experiment.HPreversal experiment.PL];
        end
        
    else
        csc=electrodes;
    end
    
    Electrodes=csc;
    
    for csc=1:length(Electrodes)
        CSC=Electrodes(csc);
        
        % load signal
        File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(CSC),'.ncs');
        FName=strcat(Path.output,filesep,'FullSignalPower',filesep,filename,filesep,'Osc','_CSC',num2str(CSC),'.mat');
        
        if exist (  FName  ) && ~flag_recalculate
            disp ('Exist, no recalculation')
        else
            mkdir(strcat(Path.output,filesep,'FullSignalPower',filesep,filename))
            cd(strcat(Path.output,filesep,'FullSignalPower',filesep,filename))
            
            if nargin < 5
                [T, Samples, samplingrate] = load_nlx(File,ExtractModeArray);
            else
                [T, Samples, samplingrate] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
            end
            
            TimeStamps=T.TimeStamps;
            [~, Samples, ~] = filter_downsample(TimeStamps, Samples, samplingrate, DownsampleFactor);
            signal=ZeroPhaseFilter(Samples,fs,[4 90]);
            
            if flag_CutData
                if length(signal)>15*60*fs
                    signal=signal(2*60*fs:15*60*fs);
                else
                    signal=signal(2*60*fs:end);
                end
            end
            
            [Power,fre,~] = pWelchSpectrum(signal,params.win_pwelch,0,2024,fs,100);
            
            save(strcat('CSC',num2str(CSC)),'Power','fre','fs');
        end
    end
end





