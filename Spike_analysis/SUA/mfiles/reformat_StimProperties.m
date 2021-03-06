function reformat_StimProperties(Experiment, StimRegion,Path)

for i=1:length(Experiment)
    experiment=Experiment(i);
    experiment.animal_ID;
    path_file=strcat(Path.temp,'\Opto_stimProperties\', experiment.name,'_Stim',StimRegion);
    
    try
        stim1=load(strcat(path_file,'_ramp','.mat'));
    catch
        stim1.StimulationProperties=[];
    end
    
    try
        stim2=load(strcat(path_file,'_squa','.mat'));
    catch
        stim2.StimulationProperties=[];
    end
    
    StimulationProperties_raw=[stim1.StimulationProperties;stim2.StimulationProperties];
    
    file2save=strcat(Path.output, filesep, 'StimulationPropertiesSUA', filesep, experiment.name);
    if ~exist(file2save)
        mkdir(file2save);
    end
    save(strcat(file2save,'\StimulationProperties_raw.mat'),'StimulationProperties_raw')
    
end







