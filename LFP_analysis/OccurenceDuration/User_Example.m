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
params.band=Band;
load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_list.mat')
main_function_get_oscillation_properties(Experiment,Path,params,Electrodes)
 