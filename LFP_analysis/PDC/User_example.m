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

flag_region_pair={'vCA1' 'PFC5';};
params.flag_region_pair=flag_region_pair;
load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_list.mat')
main_function_calculate_gPDC_xy(Experiment,Path,params,electrodes);
