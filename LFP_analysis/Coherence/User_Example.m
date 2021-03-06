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
flag_region_pair={'vCA1' 'PFC5';};
params.flag_region_pair=flag_region_pair;
load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_list.mat')
main_function_calculate_Coherence(Experiment,Path,params,Electrodes)   
