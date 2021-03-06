function main_function_CutGlue_Power(Experiment,Path,params,electrodes,flag_time_path)

fs=params.fs;
ExtractModeArray=params.ExtractModeArray;
DownsampleFactor=params.DownsampleFactor;
flag_CutData=params.flag_CutData;
lenWind=params.lenWind_pWelch;
% recalculate or not
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
                csc=experiment.LECProbereversal-11:experiment.LECProbereversal-6;
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
        FName=strcat(Path.output,filesep,'OscillationPower',filesep,filename,filesep,'Osc','_CSC',num2str(CSC),'.mat'); 
        
        if exist (  FName  ) && ~flag_recalculate
            disp ('Exist, no recalculation')
        else
            mkdir(strcat(Path.output,filesep,'OscillationPower',filesep,filename))
            cd(strcat(Path.output,filesep,'OscillationPower',filesep,filename))
            
            if nargin < 5
                [T, Samples, samplingrate] = load_nlx(File,ExtractModeArray);
            else
                [T, Samples, samplingrate] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
            end
            
            TimeStamps=T.TimeStamps;
            [~, Samples, ~] = filter_downsample(TimeStamps, Samples, samplingrate, DownsampleFactor);
            XSamples=ZeroPhaseFilter(Samples,fs,[4 90]);
            
            %  use oscillatin event
            CSC_detectOsc=CSC;
            CSC_detectOsc=experiment.LECProbereversal-8;
      
            Rdetection=load( strcat(Path.temp,'\detection_discont_events_fs_',num2str (fs),'_thr2\',filename,'\ch',num2str (CSC_detectOsc),'\oscDetect.mat') );
            StartEnd_samp=Rdetection.oscDetect; % points
            StartEnd_samp=StartEnd_samp';
            
            sil_StartEnd_samp=Rdetection.silDetect; % points
            sil_StartEnd_samp=sil_StartEnd_samp';
            
            SO{1}=StartEnd_samp;
            SO{2}=sil_StartEnd_samp;
            
            for so=1:2
                s=SO{so};
                
                % cut signal or not
                if flag_CutData
                    EndTime=s(:,2)/fs; % s
                    s=s(2:max(find(EndTime<15*60)),:);
                end
                
                % cut and glue signal
                [CutGlu_signal,~,~]=cutandglue(params,XSamples,s);
                yCutGlu=CutGlu_signal.xn;
                
                % start calculating power---------------------------
                num_segs_one_window=ceil(lenWind/params.win_sec);
                num_winds=floor(2*size(yCutGlu,1)/num_segs_one_window)-1; % sliping window, 50 overlap
                
                OscPower=[];
                for wind=1:num_winds
                    S=floor((wind-1)*num_segs_one_window/2)+1;
                    E=ceil((wind+1)*num_segs_one_window/2);
                    
                    y=[]; % glue all the segs
                    for seg=S:E
                        y=[y,yCutGlu(seg,:)];
                    end
                    
                    [oscpower,f,~] = pWelchSpectrum(y,params.win_sec,0,2024,fs,100);
                    OscPower=[OscPower;oscpower];
                end
                
                Power_SlideWind_mean=mean(OscPower,1);
                Power_SlideWind=OscPower;
                
                LFP=[];
                for seg=1:(size(yCutGlu,1))
                    y=yCutGlu(seg,:);
                    LFP=[LFP,y];
                end
                [Power_FullSig,f,~] = pWelchSpectrum(LFP,params.win_sec,0,2024,fs,100);
                
                if so==1
                    flag='Osc';
                else
                    flag='Sil';
                end
                save(strcat(flag,'_CSC',num2str(CSC)),'OscPower','f','Power_SlideWind','Power_FullSig','Power_SlideWind_mean')
                
            end
        end
        
    end
end




