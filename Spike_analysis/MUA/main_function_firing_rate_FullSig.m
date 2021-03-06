function main_function_firing_rate_FullSig(Experiment,Path,param,electrodes,flag_time_path)
% full signal , 2019.2.19, Xiaxia
flag_CutData=param.flag_CutData;
ExtractModeArray=param.ExtractModeArray;
threshold=param.th_MUA;
flag_Recalculate=param.flag_Recalculate;
FieldName=fieldnames (param);

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
            csc=experiment.LECProbereversal-12:experiment.LECProbereversal-5;
            %Electrodes=[Electrodes csc];
            Electrodes=[csc];
    end
    
    
    if ~isempty(Experiment(iExperiment).ErrorChannels)
        errorChans=intersect(Electrodes, Experiment(iExperiment).ErrorChannels);
        Electrodes(find(ismember(Electrodes,errorChans)))=[];
    end
    
    numchannels=length(Electrodes);
    
    for ichan=1:numchannels
        CSC=Electrodes(ichan);
        
        if ~exist( strcat(Path.output,filesep,'Firing_rate',filesep,filename,filesep,'CSC',num2str(CSC),'.mat')) || flag_Recalculate
            
            %load signal, time vector and sampling rate
            File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(CSC),'.ncs');
            if nargin < 5
                [~, recordingRaw,samplingrate_MUA] = load_nlx(File,ExtractModeArray);
            else
                [~, recordingRaw,samplingrate_MUA] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
            end
            
            samplingrate_MUA = round(samplingrate_MUA);
            if flag_CutData>1
                recordingRaw=recordingRaw(2*60*samplingrate_MUA:flag_CutData*60*samplingrate_MUA);
            end
            if flag_CutData==1
                recordingRaw=recordingRaw(2*60*samplingrate_MUA:14*60*samplingrate_MUA);
            end
            
            recordingMUA = ZeroPhaseFilter(recordingRaw,samplingrate_MUA,[500 5000]);
            thr = std(recordingMUA)*threshold;
            [peakLoc, ~] = peakfinderOpto(recordingMUA,-thr/2 ,-thr,-1);
            time_len=length(recordingRaw)/samplingrate_MUA;
            Firing_rate=length(peakLoc)./time_len;
            Firing_interval=diff(peakLoc)/samplingrate_MUA;
            
            % save the results
            if ~exist(  strcat(strcat(Path.output,filesep,'Firing_rate',filesep,filename))  )
                mkdir(strcat(Path.output,filesep,'Firing_rate',filesep,filename))
            end
            % mkdir(strcat(Path.output,filesep,'Firing_rate',filesep,filename))
            cd(strcat(Path.output,filesep,'Firing_rate',filesep,filename))
            save(strcat('CSC',num2str(CSC)),'Firing_rate','Firing_interval')
        end
        
    end
end






