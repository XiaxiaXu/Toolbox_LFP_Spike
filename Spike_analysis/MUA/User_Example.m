%% firing rate
clear
param.ExtractModeArray=[];param.th_MUA=5;Electrodes=1:48;param.flag_CutData=15;param.flag_Recalculate=0;
load('E:\Projects\Project-LEC\Experiment-LEC-HP-PFC_baseline_behavior\Experiment list\Experiment_list.mat')
main_function_firing_rate_FullSig(Experiment,Path,param,Electrodes);

