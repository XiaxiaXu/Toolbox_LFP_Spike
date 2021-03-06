function main_function_calculate_gPDC_xy(Experiment,Path,params,electrodes,flag_time_path)
%written by Xiaxia
%haar db sym coif bior

fs=params.fs;
th_r=params.th_r;
ExtractModeArray=params.ExtractModeArray;
DownsampleFactor=params.DownsampleFactor;
band=params.band;
flag_CutData=params.flag_CutData;
flag_mvarresidue=0; % test for the residue
wname='db4';N=4;% delete noise
nFreqs=1024;metric='diag';maxIP=50;alg=1;alpha=0.05;criterion=1;
f=fs/(2*nFreqs):fs/(2*nFreqs):fs/2;

for iExperiment=1:length(Experiment)
    iExperiment
    XX=[];
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
        if strcmp(params.flag_region_pair{1}, 'LEC2')
            chan1=experiment.LECProbereversal-11;
        end
        if strcmp(params.flag_region_pair{1}, 'LEC5')
            chan1=experiment.LECProbereversal-8;
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
        
        %-------------------------------------
        if strcmp(params.flag_region_pair{2}, 'LEC2')
            chan2=experiment.LECProbereversal-11;
        end
        if strcmp(params.flag_region_pair{2}, 'LEC5')
            chan2=experiment.LECProbereversal-8;
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
    
    simOsc_StartEnd=load( strcat(Path.Results,'\simOsc_fs_',num2str(fs),'_thr',num2str(th_r),'\',filename,'_HP_PFC.mat') );
    
%      
%     try
%         simOsc_StartEnd=load( strcat(Path.Results,'\simOsc_fs_',num2str(fs),'_thr',num2str(th_r),'\',filename,'_',num2str (chan1),'_',num2str(chan2),'.mat') );
%     catch
%         simOsc_StartEnd=load( strcat(Path.Results,'\simOsc_fs_',num2str(fs),'_thr',num2str(th_r),'\',filename,'_',num2str (chan2),'_',num2str(chan1),'.mat') );
%     end
    
    StartEnd_samp=simOsc_StartEnd.simOsc_StartEnd_samp;
    
    
    if flag_CutData
        EndTime=StartEnd_samp(:,2)/fs; % s
        StartEnd_samp=StartEnd_samp(1:max(find(EndTime<15*60)),:);
    end
    
    % load signal1----HPreversal
    CSC= chan1;
    if nargin < 5
        File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(CSC),'.ncs');
        [T, XSamples, ffs] = load_nlx(File,ExtractModeArray);
        TimeStamps=T.TimeStamps;
        [~, xSamples, fs] = filter_downsample(TimeStamps, XSamples, ffs, DownsampleFactor);
    else
        [T, XSamples, ffs] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
        TimeStamps=T.TimeStamps;
        [~, xSamples, fs] = filter_downsample(TimeStamps, XSamples, ffs, DownsampleFactor);
    end
    xSamples=ZeroPhaseFilter(xSamples,fs,[band(1) band(2)]);
    [CutGlu_X,~,~]=cutandglue(params,xSamples,StartEnd_samp);
    XCutGlu=CutGlu_X.xn;
    
    % load signal2 ----------PFC
    
    num=size(chan2,1);
    
    for chan=1:num
        CSC=chan2(chan);
        if nargin < 5
            File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(CSC),'.ncs');
            [T, YSamples, ffs] = load_nlx(File,ExtractModeArray);
            TimeStamps=T.TimeStamps;
            [~, ySamples, fs] = filter_downsample(TimeStamps, YSamples, ffs, DownsampleFactor);
        else
            [T, YSamples, ffs] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
            TimeStamps=T.TimeStamps;
            [~, ySamples, fs] = filter_downsample(TimeStamps, YSamples, ffs, DownsampleFactor);
        end
        
        ySamples=ZeroPhaseFilter(ySamples,fs,[band(1) band(2)]);
        [CutGlu_Y,~,~]=cutandglue(params,ySamples,StartEnd_samp);
        YCutGlu=CutGlu_Y.xn;
        
        % calculating gPDC
        
        c12=[];c21=[];
        for seg=1:(size(YCutGlu,1))
            x=XCutGlu(seg,:);
            y=YCutGlu(seg,:);
            
            x=detrend(x);
            y=detrend(y);
            
            [THR,SORH,KEEPAPP] = ddencmp('den','wp',x);
            x=wdencmp('gbl',x, wname,N,THR,SORH,KEEPAPP) ;
            
            [THR,SORH,KEEPAPP] = ddencmp('den','wp',y);
            y=wdencmp('gbl',y, wname,N,THR,SORH,KEEPAPP) ;
            
            u=[x' y'];
            c=PDC_computation(u,nFreqs,metric,maxIP,alg,alpha,criterion,flag_mvarresidue);
            
            if ~isempty(c)
                c12=[c12,c.c12];
                c21=[c21,c.c21];
            end
        end
        c12(find(isnan(c12)))=0;
        c21(find(isnan(c21)))=0;
        c12_mean=mean(c12,2);
        c21_mean=mean(c21,2);
        
        if ~exist(  (strcat(Path.output,filesep,'gPDC',filesep,filename))  )
            mkdir( (strcat(Path.output,filesep,'gPDC',filesep,filename))  )
        end
        
        cha1=chan1;
        cha2=chan2(chan);
        cd( (strcat(Path.output,filesep,'gPDC',filesep,filename)) )
        save( strcat('XLFP',num2str(cha1),'_YLFP',num2str(cha2)),'c12' ,'c21','c12_mean','c21_mean','fs','f');
        
    end
end




