function main_function_OscViewer_ContinuousData(Experiment,iExperiment, params,electrodes,StartWindow, flag_time_path)
DownsampleFactor=params.DownsampleFactor;
ExtractModeArray=params.ExtractModeArray;
window_length=params.window_length;
flag_CutData=params.flag_CutData;
mn=3;

experiment=Experiment(iExperiment);
if electrodes==1
    Electrodes=experiment.HPreversal ;
elseif electrodes==2
    Electrodes=experiment.PL;
else
    Electrodes=electrodes;
end

CSC=Electrodes;
%load signal, time vector and sampling rate
File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(CSC),'.ncs');
if nargin < 6
    [T, Samples, Fs] = load_nlx(File,ExtractModeArray);
else
    [T, Samples, Fs] = load_nlx_stimulation_baseline(experiment, CSC,flag_time_path);
end
Fs=round(Fs);
TimeStamps=T.TimeStamps;
[~, samples, fs] = filter_downsample(TimeStamps, Samples,Fs, DownsampleFactor);
signal=ZeroPhaseFilter(samples,fs,[1 100]);
signal_4_100=ZeroPhaseFilter(samples,fs,[4 100]);
signal_4_30=ZeroPhaseFilter(samples,fs,[4 30]);
signal_30_100=ZeroPhaseFilter(samples,fs,[30 100]);
signal_MUA=ZeroPhaseFilter(Samples,Fs,[500 1000]);

if flag_CutData
    if length(signal)>20*60*fs
        signal=signal(2*60*fs:20*60*fs);
        signal_4_100=signal_4_100(2*60*fs:20*60*fs);
        signal_4_30=signal_4_30(2*60*fs:20*60*fs);
        signal_30_100=signal_30_100(2*60*fs:20*60*fs);
        sig_MUA=signal_MUA(2*60*Fs:20*60*Fs);
    else
        signal=signal(2*60*fs:end);
        signal_4_100=signal_4_100(2*60*fs:end);
        signal_4_30=signal_4_30(2*60*fs:end);
        signal_30_100=signal_30_100(2*60*fs:end);
        sig_MUA=signal_MUA(2*60*Fs:end);
    end
end

%Loop through each window
Sig=signal;
lenWind=window_length*fs;
num_winds=floor(2*length(Sig)/lenWind)-5; % sliping window, 50 overlap

jj = StartWindow: num_winds;
f = figure('Color', 'white');
hold on
colormap('jet')
for num=jj(2:end)
    
    lenWind=window_length*fs;
    S=(num-1)*lenWind/2+1;
    E=(num+1)*lenWind/2;
    
    disp(['Window: ' num2str(num) '/' num2str(num_winds)]);
    figure(f);
    
    ii=0;
%     ii=ii+1;subplot(mn,1,ii); 
%     LFP=signal(S:E);plot(LFP);set(gca,'xtick',[]); %ylim([-100 100]);
    
    ii=ii+1;subplot(mn,1,ii);
    LFP=signal_4_100(S:E);plot(LFP);%set(gca,'xtick',[]);%ylim([-100 100]);
    
%     ii=ii+1;subplot(mn,1,ii);
%     LFP=signal_30_100(S:E);plot(LFP);set(gca,'xtick',[]);%ylim([-100 100]);
    
    ii=ii+1;subplot(mn,1,ii);
    LFP=signal(S:E);time=1/fs:1/fs:window_length; % ms
    [waveLFP, periodLFP, ~, ~] = wt([time; LFP],'S0',1/100,'MaxScale',1/2);
    powerLFP      = (abs(waveLFP)).^2 ;
    sigma2LFP     = var(LFP);
    imagesc(time,log2(periodLFP),log2(abs(powerLFP/sigma2LFP)));
%     contourf(time,log2(periodLFP),log2(abs(powerLFP/sigma2LFP)));
% image(time,log2(periodLFP),log2(abs(powerLFP/sigma2LFP)));

    %colorbar
    clim=get(gca,'clim'); %center color limits around log2(1)=0
    globclim4_100=max(clim(2),3); %global
    clim=[0 1]*globclim4_100;
    set(gca,'clim',clim)
    Yticks = 2.^(fix(log2(min(periodLFP))):fix(log2(max(periodLFP))));
    set(gca,'YLim',log2([min(periodLFP),max(periodLFP)]), ...
        'YDir','reverse', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(1./Yticks'), ...
        'layer','top')
    title('1-100Hz Wavelet')
    set(gca,'xtick',[]) % Removes X axis label
    
    % for MUA
   ii=ii+1;subplot(mn,1,ii);
    lenWind=window_length*Fs;
    S=(num-1)*lenWind/2+1;
    E=(num+1)*lenWind/2;
    MUA=sig_MUA(S:E);
    plot(MUA);
    %ylim([-35 35])
    
    try
        reply1=input('Next event (Enter)');
    catch
        reply1=input('Next event (Enter)');
    end
end

