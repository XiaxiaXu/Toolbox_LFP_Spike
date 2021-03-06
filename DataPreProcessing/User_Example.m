%% detect Events
clear
downsampleFactor=[32]; % fs=1000;
DownsampleFactor=downsampleFactor(1);    
band=[4 30]; % theta & beta
param=  struct('tapers',[5 9], 'pad', 0,'fs', 32000/DownsampleFactor,'fpass',band,'err',[1 0.05],'trialave',0);
param.win_sec=1; param.th_r=2;param.DownsampleFactor=DownsampleFactor;param.ExtractModeArray=[];
Electrodes=[];
param.flag_region={'LEC'};
param.flag_cutSignal=15;

load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_list.mat')
main_function_detectEvents_new(Experiment,Path,band,param,Electrodes);

param.flag_region_pair= {'vCA1' 'PFC5'};
main_function_getSimEvents_new(Experiment,Path,param,Electrodes);
