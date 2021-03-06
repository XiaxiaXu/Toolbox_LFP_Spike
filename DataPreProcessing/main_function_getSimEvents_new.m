function  main_function_getSimEvents_new(Experiment,Path,param,electrodes)
%written by Sabine Gretenkord - last modified 25.10.2016
%finds logical vector and time stamps for simultaneous oscillatory periods
%also gives out logical vector and time stamps for non-oscillatory periods
%as well as logical vector for oscillatory periods occurring in region 1 but not region
%2 and vice versa
th_r=param.th_r;
fs=param.fs;
num_sd=param.th_r;
for iExperiment=1:length(Experiment)
    iExperiment
    filename=Experiment(iExperiment).name;
    resultFile=strcat(Path.output,filesep,'simOsc_fs_',num2str(fs),'_thr',num2str(th_r),filesep,filename,'.mat');
    experiment=Experiment(iExperiment);
    FieldName=fieldnames (param);
    if ismember('flag_region_pair',FieldName)
        if strcmp(param.flag_region_pair{1}, 'LEC2')
            Electrodes(1)=experiment.LECProbereversal-7;
        end
        
        if strcmp(param.flag_region_pair{1}, 'LEC5')
            Electrodes(1)=experiment.LECProbereversal-4;
        end
        
        if strcmp(param.flag_region_pair{1}, 'vCA1')
            Electrodes(1)=experiment.HPreversal;
        end
        if strcmp(param.flag_region_pair{1}, 'PFC2')
            Electrodes(1)=Experiment(iExperiment).PL;
        end
         if strcmp(param.flag_region_pair{1}, 'PFC5')
            Electrodes(1)=Experiment(iExperiment).PL+12;
        end
        %--------------------------------------
        if strcmp(param.flag_region_pair{2}, 'LEC5')
            Electrodes(2)=experiment.LECProbereversal-4;
        end
        
        if strcmp(param.flag_region_pair{2}, 'LEC2')
            Electrodes(2)=experiment.LECProbereversal-7;
        end
        
        if strcmp(param.flag_region_pair{2}, 'vCA1')
            Electrodes(2)=experiment.HPreversal;
        end
        if strcmp(param.flag_region_pair{2}, 'PFC2')
            Electrodes(2)=Experiment(iExperiment).PL;
        end
        if strcmp(param.flag_region_pair{2}, 'PFC5')
            Electrodes(2)=Experiment(iExperiment).PL+12;
        end
        
    elseif isempty(electrodes)
        Electrodes=[Experiment(iExperiment).HPreversal Experiment(iExperiment).PL];
    else
        Electrodes=electrodes;
    end
    
    Electrodes1=Electrodes(1);
    Electrodes2=Electrodes(2);
    
    if exist(resultFile)==2
        
        disp( strcat( 'Experiment ',mat2str(iExperiment),': "',Experiment(iExperiment).name,'" sim Osc exists')  )
        
        %else
    end
    
    %load event detection region 1
    
    eventDetectionFolder1=strcat(Path.temp,filesep,'detection_discont_events_fs_',num2str(fs),'_thr',num2str(num_sd),'\',filename,'\');
    cd(strcat(eventDetectionFolder1,filesep))
    Region1=load(strcat(eventDetectionFolder1,'ch',num2str(Electrodes1),'\oscDetect'));
    
    %load event detection region 2
    eventDetectionFolder2=strcat(Path.temp,filesep,'detection_discont_events_fs_',num2str(fs),'_thr',num2str(num_sd),'\',filename,'\');
    cd(strcat(eventDetectionFolder2,filesep))
    Region2=load(strcat('ch',num2str(Electrodes2),'\oscDetect'));
    
    %one-region oscillations
    Region1OnlyOscLogi= Region1.oscLogi &  ~Region2.oscLogi; %logical
    Region2OnlyOscLogi= ~Region1.oscLogi &  Region2.oscLogi; %logical
    
    %two-region oscillations: timestamps (samples)
    simOscLogi= Region1.oscLogi &  Region2.oscLogi; %logical
    transLogi=[diff(simOscLogi),0];
    simOscStart_samp=find(transLogi==1);
    simOscEnd_samp=find(transLogi==-1);
    simOsc_StartEnd_samp=[simOscStart_samp',simOscEnd_samp'];
    
    %oscillations in neither of the two regions: timestamps (samples)
    noOscLogi= ~Region1.oscLogi &  ~Region2.oscLogi;%logical
    transLogi=[diff(noOscLogi),0];
    noOscStart_samp=find(transLogi==1);
    noOscEnd_samp=find(transLogi==-1);
    noOsc_StartEnd_samp=[noOscStart_samp(1:end-1)',noOscEnd_samp(2:end)'];
    
    %load time vector (of LFP)
    %chanChoice=Electrodes2;
    %File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(chanChoice),'.ncs');
    %[TimeStamps, Samples, fs] = load_nlx(File,ExtractModeArray);
    %[TimeStamps, Samples, fs] = filter_downsample(TimeStamps, Samples, fs, DownsampleFactor);
    %T=TimeStamps/10^3; %time in seconds
    
    T=1/fs:1/fs:10000;
    
    %two-region oscillations: calculate timestamps (seconds)
    simOscStart_time=T(simOscStart_samp);
    simOscEnd_time=T(simOscEnd_samp);
    simOsc_StartEnd_time=[simOscStart_time',simOscEnd_time'];
    
    %oscillations in neither of the two regions: timestamps (seconds)
    noOscStart_time=T(noOscStart_samp);
    noOscEnd_time=T(noOscEnd_samp);
    noOsc_StartEnd_time=[noOscStart_time',noOscEnd_time'];
    
    %save  result and figure
    mkdir(strcat(Path.output,filesep,'simOsc_fs_',num2str(fs),'_thr',num2str(th_r)))
    cd(strcat(Path.output,filesep,'simOsc_fs_',num2str(fs),'_thr',num2str(th_r)))
%     save(filename,'simOsc_StartEnd_time','simOsc_StartEnd_samp', 'noOsc_StartEnd_time','noOsc_StartEnd_samp','simOscLogi','noOscLogi','Region1OnlyOscLogi','Region2OnlyOscLogi')
    save(strcat(filename,'_',num2str(Electrodes1),'_',num2str(Electrodes2)),'simOsc_StartEnd_time','simOsc_StartEnd_samp', 'noOsc_StartEnd_time','noOsc_StartEnd_samp','simOscLogi','noOscLogi','Region1OnlyOscLogi','Region2OnlyOscLogi')
   
end

clearvars -except Experiment iExperiment Path th_r  param  Electrodes DownsampleFactor ExtractModeArray fs

end




