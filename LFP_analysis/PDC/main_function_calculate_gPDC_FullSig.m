function main_function_calculate_gPDC_FullSig(Experiment,Path,params,electrodes,flag_time_path)
%haar db sym coif bior

DownsampleFactor=params.DownsampleFactor;
fs=params.fs;
ExtractModeArray=params.ExtractModeArray;
window_length=params.window_length;
flag_CutData=params.flag_CutData;

wname='db4';N=4;% reduce noise
nFreqs=512;metric='diag';maxIP=50;alg=1;criterion=1;
%alpha=0.05; 
alpha=0; 
f=fs/(2*nFreqs):fs/(2*nFreqs):fs/2;

for iExperiment=1:length(Experiment)
    experiment=Experiment(iExperiment);
    iExperiment
    filename=Experiment(iExperiment).name;
    if isempty(electrodes)
        Electrodes=[experiment.HPreversal experiment.PL];
    else
        Electrodes=electrodes;
    end
    
    numchannels=length(Electrodes);
    Sig=[];
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
        Fs=round(Fs);
        TimeStamps=T.TimeStamps;
        [~, samples, fs] = filter_downsample(TimeStamps, Samples,Fs, DownsampleFactor);
        signal=ZeroPhaseFilter(samples,fs,[3 45]);
        
        if flag_CutData
            if length(signal)>20*60*fs
                signal=signal(2*60*fs:20*60*fs);
            else
                signal=signal(2*60*fs:end);
            end
        end
        
        Sig(ichan,:)=signal;
    end
    
    lenWind=window_length*fs;
    num_winds=floor(2*length(Sig(1,:))/lenWind)-5; % sliping window, 50 overlap
    % start calculating
    c12=[];c21=[];
    for num=1:num_winds
        S=(num-1)*lenWind/2+1;
        E=(num+1)*lenWind/2;
        x=Sig(1,S:E);
        y=Sig(2,S:E);
        x=detrend(x);
        y=detrend(y);
        
        [THR,SORH,KEEPAPP] = ddencmp('den','wp',x);
        x=wdencmp('gbl',x, wname,N,THR,SORH,KEEPAPP) ;
        
        [THR,SORH,KEEPAPP] = ddencmp('den','wp',y);
        y=wdencmp('gbl',y, wname,N,THR,SORH,KEEPAPP) ;
        
        u=[x' y'];
        flag_mvarresidue=1;
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
    
    cd( (strcat(Path.output,filesep,'gPDC',filesep,filename)) )
    save( strcat('XLFP',num2str(Electrodes(1)),'_YLFP',num2str(Electrodes(2))),'c12' ,'c21','c12_mean','c21_mean','fs','f','window_length');
end


