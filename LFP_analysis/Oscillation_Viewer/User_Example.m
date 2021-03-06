clc
clear
DownsampleFactor=32;
params.DownsampleFactor=DownsampleFactor;
params.ExtractModeArray=[];
params.flag_CutData=1; % 
params.window_length=100;
BrainRegion=2; % 1--HP; 2--PFC; else real number
StartWindow=2; 
load('D:\Project-DISC1 in PFC-Neonatal\Experiment list_Analysis\Experiment_Henrike_Con.mat')
load('D:\Project-DISC1 in PFC-Neonatal\Experiment list_Analysis\Experiment_Xiaxia_Con.mat')

iExperiment=12;
main_function_OscViewer_ContinuousData(Experiment,iExperiment, params,BrainRegion,StartWindow)