function [] = main_function_downsample_OptoStim(Experiment,Path)

%written by Sabine Gretenkord, last modified 21.08.2017
%downsample the stimulation signal to align it like the LFP
%to make sure the segments are chosen correctly

for iExperiment=1:length(Experiment)
    
    %get Opto Stim Data data------------
    folder_stim=strcat(Path.temp,'\nlx_load_OptoStim\',Experiment(iExperiment).name );
    cd(folder_stim)
    
    if ~exist('STIM1A_downsampled.mat')
        
        Stim=load( 'STIM1A.mat') ;
        
        signal=Stim.signal(1:20:end);
        time=Stim.time(1:20:end);
        fs=Stim.fsInput/20;
        
        save('STIM1A_downsampled','signal','time','fs')
        
    end

    %     figure
    %     hold on
    %     plot(Stim.time,Stim.signal,'k')
    %     plot(time,signal-1000,'r')
    
end

end




