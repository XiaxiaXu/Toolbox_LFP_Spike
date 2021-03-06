%% Expect running time is less than 5 min.

%% detect Events
clear
downsampleFactor=[32 160]; % fs=1000;
for i=1:length(downsampleFactor)
    DownsampleFactor=downsampleFactor(i);
    
    band=[4 30]; % theta & beta
    param=  struct('tapers',[5 9], 'pad', 0,'fs', 32000/DownsampleFactor,'fpass',band,'err',[1 0.05],'trialave',0);
    param.win_sec=1; param.th_r=2;param.DownsampleFactor=DownsampleFactor;param.ExtractModeArray=[];
    Electrodes=[];
    param.flag_region={'LEC'};
    param.flag_cutSignal=15;
    param.flag_recalculate=0;
    
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
    main_function_detectEvents_new(Experiment,Path,band,param,Electrodes);
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
    Experiment(33:34)=[];
    main_function_detectEvents_new(Experiment,Path,band,param,Electrodes);
    
    flag_region_pair={'LEC2','vCA1';'LEC5','vCA1';...
        'LEC5','PFC2';'LEC5','PFC5';'LEC2','PFC2';'LEC2','PFC5';...
        'vCA1' 'PFC2';'vCA1' 'PFC5';};
    
    for j=1:size (flag_region_pair,1)
        j
        param.flag_region_pair= flag_region_pair(j,:);
        load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
        main_function_getSimEvents_new(Experiment,Path,param,Electrodes);
        load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
        main_function_getSimEvents_new(Experiment,Path,param,Electrodes);
    end
    
end

%% Oscillation properties- Occurance+Amplitude+Duration----------HP and PFC
clear

params.DownsampleFactor=32;
params.th_LFP=2; % 2~2.5
params.ExtractModeArray=[];
params.flag_CutData=1;
Electrodes=[];
Band=[4 30];
params.flag_band=[];
params.flag_LEC='Yes';

for i=1:size(Band,1)
    params.band=Band(i,:);
    
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
    main_function_get_oscillation_properties(Experiment,Path,params,Electrodes)
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
    main_function_get_oscillation_properties(Experiment,Path,params,Electrodes)
end

%% Oscillation Power
clear
DownsampleFactor=32; % fs=1000;
param.fs=32000/DownsampleFactor;
param.win_sec=1; param.th_r=2;param.DownsampleFactor=DownsampleFactor;
param.ExtractModeArray=[];param.flag_CutData=1;
param.lenWind_pWelch=10;
flag_region={'LEC2','LEC5',...
    'PFC2','PFC5',...
    'vCA1'};

param.flag_recalculate=0;
Electrodes=[];
load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
for j=1:length (flag_region)
    param.flag_region=flag_region{j};
    main_function_CutGlue_Power(Experiment,Path,param,Electrodes)
end
Electrodes=[];
load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
for j=1:length (flag_region)
    param.flag_region=flag_region{j};
    main_function_CutGlue_Power(Experiment,Path,param,Electrodes)
end

%% Coherence
clear
DownsampleFactor=32;
params=  struct('tapers',[5 9], 'pad', 0,'Fs', 32000/DownsampleFactor,'fpass',[1 100],'err',[1 0.05],'trialave',0);
params.DownsampleFactor=DownsampleFactor;
params.fs=1000;
params.th_LFP=2;
params.ExtractModeArray=[];
params.flag_CutData=1; %
params.band=[4 90];
params.win_sec=1;
params.lenWind_pWelch=10;
params.flag_recalculate=1;
Electrodes=[];
flag_region_pair={'LEC2','vCA1';'LEC5','vCA1';...
    'LEC5','PFC2';'LEC5','PFC5';'LEC2','PFC2';'LEC2','PFC5';...
    'vCA1' 'PFC2';'vCA1' 'PFC5';};

flag_region_pair={'LEC2','vCA1';'LEC2','PFC5';'vCA1' 'PFC5';};

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
    main_function_calculate_Coherence(Experiment,Path,params,Electrodes)
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
    main_function_calculate_Coherence(Experiment,Path,params,Electrodes)
end

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
    main_function_calculate_Coherence_Surrogate(Experiment,Path,params,Electrodes)

    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
    main_function_calculate_Coherence_Surrogate(Experiment,Path,params,Electrodes)
end

%% Coherence_ref
clear
DownsampleFactor=32;
params=  struct('tapers',[5 9], 'pad', 0,'Fs', 32000/DownsampleFactor,'fpass',[1 100],'err',[1 0.05],'trialave',0);
params.DownsampleFactor=32;
params.fs=1000;
params.th_LFP=2;
params.ExtractModeArray=[];
params.flag_CutData=1; %
params.band=[4 90];
params.win_sec=1;
params.lenWind_pWelch=10;
params.flag_recalculate=0;
Electrodes=[];
flag_region_pair={'LEC2','vCA1';'LEC5','vCA1';...
    'LEC5','PFC2';'LEC5','PFC5';'LEC2','PFC2';'LEC2','PFC5';...
    'vCA1' 'PFC2';'vCA1' 'PFC5';};

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_Con.mat')
    Experiment(6)=[];
    main_function_calculate_Coherence_Re_reference(Experiment,Path,params,Electrodes)
    
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_GE.mat')
    Experiment([17 19])=[];
    main_function_calculate_Coherence_Re_reference(Experiment,Path,params,Electrodes)
end

%% cross correlation
clear;
DownsampleFactor=32;
params=  struct('tapers',[5 9], 'pad', 0,'Fs', 32000/DownsampleFactor,'fpass',[1 100],'err',[1 0.05],'trialave',0);
params.DownsampleFactor=32;
params.fs=32000/DownsampleFactor;
params.th_LFP=2;
params.ExtractModeArray=[];
params.flag_CutData=1; %
params.win_sec=1;
params.Band=[4 12;12 30;]; % it is unreasonable to calculate the directionality at gamma rhythm
params.flag_recalculate=0;
electrodes=[];
flag_region_pair={'LEC2','vCA1';'LEC5','vCA1';...
    'LEC5','PFC2';'LEC5','PFC5';'LEC2','PFC2';'LEC2','PFC5';...
    'vCA1' 'PFC2';'vCA1' 'PFC5';};

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_Con.mat')
    main_function_calculate_Correlation_CutGlue(Experiment,Path,params,electrodes)
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_GE.mat')
    main_function_calculate_Correlation_CutGlue(Experiment,Path,params,electrodes)
end

%% cross correlation_Ref
clear;
DownsampleFactor=32;
params=  struct('tapers',[5 9], 'pad', 0,'Fs', 32000/DownsampleFactor,'fpass',[1 100],'err',[1 0.05],'trialave',0);
params.DownsampleFactor=32;
params.fs=32000/DownsampleFactor;
params.th_LFP=2;
params.ExtractModeArray=[];
params.flag_CutData=1; %
params.win_sec=1;
params.Band=[4 12;12 30;]; % it is unreasonable to calculate the directionality at gamma rhythm
params.flag_recalculate=0;
electrodes=[];
flag_region_pair={'LEC2','vCA1';'LEC5','vCA1';...
    'LEC5','PFC2';'LEC5','PFC5';'LEC2','PFC2';'LEC2','PFC5';...
    'vCA1' 'PFC2';'vCA1' 'PFC5';};

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_Con.mat')
    Experiment(6)=[];
    main_function_calculate_Correlation_CutGlue_Re_reference(Experiment,Path,params,electrodes)
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_GE.mat')
    Experiment([17 19])=[];
    main_function_calculate_Correlation_CutGlue_Re_reference(Experiment,Path,params,electrodes)
end

%% gPDC
clear
DownsampleFactor=160;
params.fs=32000/DownsampleFactor;
params.DownsampleFactor=DownsampleFactor;
params.ExtractModeArray=[];
params.th_r=2;
params.band=[4 45];
params.win_sec=1;
params.flag_CutData=1;
electrodes=[];

flag_region_pair={'LEC2','vCA1';'LEC5','vCA1';...
    'LEC5','PFC2';'LEC5','PFC5';'LEC2','PFC2';'LEC2','PFC5';...
    'vCA1' 'PFC2';'vCA1' 'PFC5';};

flag_region_pair={'LEC','PFC'; 'vCA1','PFC';'LEC','vCA1' };

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
    main_function_calculate_gPDC_xy(Experiment,Path,params,electrodes);
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
    main_function_calculate_gPDC_xy(Experiment,Path,params,electrodes);
end


%% gPDC_reference
clear
DownsampleFactor=160;
params.fs=32000/DownsampleFactor;
params.DownsampleFactor=DownsampleFactor;
params.ExtractModeArray=[];
params.th_r=2;
params.band=[4 45];
params.win_sec=1;
params.flag_CutData=1;
electrodes=[];

flag_region_pair={'LEC2','vCA1';'LEC5','vCA1';...
    'LEC5','PFC2';'LEC5','PFC5';'LEC2','PFC2';'LEC2','PFC5';...
    'vCA1' 'PFC2';'vCA1' 'PFC5';};
reference_chan=35; % LEC channel

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_Con.mat')
    Experiment(6)=[];
    main_function_calculate_gPDC_xy_Re_reference(Experiment,Path,params,electrodes,reference_chan);
    load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_Con.mat')
    Experiment([17 19])=[];
    main_function_calculate_gPDC_xy_Re_reference(Experiment,Path,params,electrodes,reference_chan);
end

%% firing rate
clear
param.ExtractModeArray=[];param.th_MUA=5;Electrodes=1:48;param.flag_CutData=15;param.flag_Recalculate=0;

load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
main_function_firing_rate_FullSig(Experiment,Path,param,Electrodes);

load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
main_function_firing_rate_FullSig(Experiment,Path,param,Electrodes);

%% CSD
clear
param.ExtractModeArray=[];param.flag_Recalculate=1; param.DownsampleFactor=32;param.band=[1 90];param.spacing=10^-6;

param.StartTime=5; % min
param.flag_CutData=60; % ms

load('D:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
Experiment=Experiment(2);
main_function_CSD(Experiment,Path,param)


%% HP sharp wave -- LEC spike
clear
clc

params.DownsampleFactor=32;
params.threshold_SharpWave=5;
params.ExtractModeArray=[];
params.StartTime_min=10;
params.EndTime_min=15;
params.flag_CutData=0;
params.th_MUA=4.5;
params.pre_post_SWR=3;
params.SpikeRegion='LEC';
params.SpikeRegion='LEC5';
params.SpikeRegion='LEC2';

load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
Experiment=Experiment(11:20);
main_function_SharpWave_Spike(Experiment,Path,params)

load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
main_function_SharpWave_Spike(Experiment,Path,params)


%% HP sharp wave -- LEC spike
clear
clc

params.DownsampleFactor=32;
params.threshold_SharpWave=5;
params.ExtractModeArray=[];
params.StartTime_min=10;
params.EndTime_min=15;
params.flag_CutData=1;

params.pre_post_SWR=3;
params.LFPRegion='LEC';
params.LFPRegion='PFC';

load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
Experiment(23)=[];
main_function_SharpWave_evokeLFP(Experiment,Path,params)

load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
main_function_SharpWave_evokeLFP(Experiment,Path,params)




%
%% SIC
clear
DownsampleFactor=32;
params.fs=32000/DownsampleFactor;
params.DownsampleFactor=DownsampleFactor;
params.ExtractModeArray=[];
params.th_r=2;
params.win_sec=1;
params.flag_CutData=1;
electrodes=[];


flag_region_pair={'LEC','PFC'; 'vCA1','PFC';'LEC','vCA1' };
%flag_region_pair={'LEC','vCA1'};

for j=1:size (flag_region_pair,1)
    params.flag_region_pair=flag_region_pair(j,:);
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_Con.mat')
    main_function_calculate_SIC(Experiment,Path,params,electrodes);
    load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_LEC-HP-PFC_All_GE.mat')
    main_function_calculate_SIC(Experiment,Path,params,electrodes);

end














