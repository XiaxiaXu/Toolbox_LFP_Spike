function main_function_calculate_Coherence(Experiment,Path,params,electrodes,flag_time_path)

fs=params.fs;
th_r=params.th_LFP;
ExtractModeArray=params.ExtractModeArray;
flag_CutData=params.flag_CutData;
DownsampleFactor=params.DownsampleFactor;
band=params.band;
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
    experiment=Experiment(iExperiment);
    filename=Experiment(iExperiment).name;
    %% prepare the channels
    if isempty(electrodes)
        Electrodes=experiment.PL;
    else
        Electrodes=electrodes;
    end
    
    FieldName=fieldnames (params);
    if ismember('flag_region_pair',FieldName)
        if strcmp(params.flag_region_pair{1}, 'LEC5')
            chan1=experiment.LECProbereversal-4;
        end
        if strcmp(params.flag_region_pair{1}, 'LEC2')
            chan1=experiment.LECProbereversal-7;
        end
        if strcmp(params.flag_region_pair{1}, 'vCA1')
            chan1=Experiment(iExperiment).HPreversal;
        end
        if strcmp(params.flag_region_pair{1}, 'PFC2')
            chan1=Electrodes;
        end
        if strcmp(params.flag_region_pair{1}, 'PFC5')
            chan1=experiment.PL+12;
        end
        %--------------------------------------------
        if strcmp(params.flag_region_pair{2}, 'LEC5')
            chan2=experiment.LECProbereversal-4;
        end
        if strcmp(params.flag_region_pair{2}, 'LEC2')
            chan2=experiment.LECProbereversal-7;
        end
        if strcmp(params.flag_region_pair{2}, 'vCA1')
            chan2=Experiment(iExperiment).HPreversal;
        end
        if strcmp(params.flag_region_pair{2}, 'PFC2')
            chan2=Electrodes;
        end
        if strcmp(params.flag_region_pair{2}, 'PFC5')
            chan2=experiment.PL+12;
        end
    else
        chan1=Experiment(iExperiment).HPreversal;
        chan2=Electrodes;
    end
    
    FName=strcat(Path.output,filesep,'Coherence',filesep,filename,filesep,'LFP',num2str(chan1),'_LFP',num2str(chan2),'.mat');  
    % save the results
    if exist(  FName )&& ~flag_recalculate
        disp ('Exist, no recalculation')
    else
        mkdir(strcat(Path.output,filesep,'Coherence',filesep,filename))
        cd(strcat(Path.output,filesep,'Coherence',filesep,filename))
        % use co-occuring LFP
        %     simOsc_StartEnd=load( strcat(Path.Results,'\simOsc_fs_',num2str(fs),'_thr',num2str(th_r),'\',filename,'.mat') );
        %     StartEnd_samp=simOsc_StartEnd.simOsc_StartEnd_samp;
        
        try
            simOsc_StartEnd=load( strcat(Path.Results,'\simOsc_fs_',num2str(fs),'_thr',num2str(th_r),'\',filename,'_',num2str (chan1),'_',num2str(chan2),'.mat') );
        catch
            simOsc_StartEnd=load( strcat(Path.Results,'\simOsc_fs_',num2str(fs),'_thr',num2str(th_r),'\',filename,'_',num2str (chan2),'_',num2str(chan1),'.mat') );
        end
        
        StartEnd_samp=simOsc_StartEnd.simOsc_StartEnd_samp;
        
        if flag_CutData==1
            EndTime=StartEnd_samp(:,2)/fs; % s
            StartEnd_samp=StartEnd_samp(1:max(find(EndTime<15*60)),:);
        elseif flag_CutData==2
            EndTime=StartEnd_samp(:,2)/fs; % s
            StartTime=StartEnd_samp(:,1)/fs; % s
            StartEnd_samp=StartEnd_samp(min(find(StartTime>EndTime-15*60)):end,:);
        end
        
        %% load signal-1
        
        File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(chan1),'.ncs');
        if nargin < 5
            [T, xSamples,samplingrate_MUA] = load_nlx(File,ExtractModeArray);
        else
            [T, xSamples,samplingrate_MUA] = load_nlx_stimulation_baseline(experiment, chan1,flag_time_path);
        end
        TimeStamps=T.TimeStamps;
        [~, XSamples, fs] = filter_downsample(TimeStamps, xSamples, samplingrate_MUA, DownsampleFactor);
        XSamples=ZeroPhaseFilter(XSamples,fs,[band(1) band(2)]);
        [CutGlu_X,~,~]=cutandglue(params,XSamples,StartEnd_samp);
        XCutGlu=CutGlu_X.xn;
        
        % load signal2 ----------
        File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(chan2),'.ncs');
        if nargin < 5
            [T, ySamples,samplingrate_MUA] = load_nlx(File,ExtractModeArray);
        else
            [T, ySamples,samplingrate_MUA] = load_nlx_stimulation_baseline(experiment, chan2,flag_time_path);
        end
        TimeStamps=T.TimeStamps;
        [~, YSamples, fs] = filter_downsample(TimeStamps, ySamples, samplingrate_MUA, DownsampleFactor);
        YSamples=ZeroPhaseFilter(YSamples,fs,[band(1) band(2)]);
        [CutGlu_Y,~,~]=cutandglue(params,YSamples,StartEnd_samp);
        YCutGlu=CutGlu_Y.xn;
        
        % ----- multi-taper------
        Coh_abs_RealImag=[];Coh_abs_Imag=[];
        for seg=1:(size(YCutGlu,1))
            x=XCutGlu(seg,:);
            y=YCutGlu(seg,:);
            [Coh,~,~,~,~,f]=coherencyc_addImag(x,y,params);
            Coh_abs_RealImag(:,seg)=Coh.abs_RealImag';
            Coh_abs_Imag(:,seg)=Coh.abs_Imag';
        end
        
        % ---- pWelch---- msCoherence
        size(YCutGlu,1)
        num_segs_one_window=ceil(lenWind/params.win_sec);
        num_winds=floor(2*size(YCutGlu,1)/num_segs_one_window)-1 ;% sliping window, 50 overlap
        
        Coh_mscohere_abs_Imag=[];
        for wind=1:num_winds
            S=floor((wind-1)*num_segs_one_window/2)+1;
            E=ceil((wind+1)*num_segs_one_window/2);
            
            x=[];y=[]; % glue all the segs
            for seg=S:E
                y=[y,YCutGlu(seg,:)];
                x=[x,XCutGlu(seg,:)];
            end
            
            [Coh_mscoherence,f_mscohere]=ImCohere(x,y,params.win_sec,0,fs,fs);
            Coh_mscohere_abs_Imag(:,wind)=Coh_mscoherence.abs_Imag;
            
        end
        
        Coh_mscohere_abs_Imag_mean=mean(Coh_mscohere_abs_Imag,2);
        save(strcat('LFP',num2str(chan1),'_LFP',num2str(chan2)),'f','Coh_abs_RealImag','Coh_abs_Imag', 'Coh_mscohere_abs_Imag_mean','Coh_mscohere_abs_Imag','f_mscohere')
    end
end


