function SUAgetStimulationProperties_adjustformat(experiments, save_data, repeatCalc,Path,Path_stimOld, StimRegion)


for n_animal = 1:length(experiments)
    n_animal
    experiment = experiments(n_animal);
    if ~isempty(experiment.animal_ID)
        if  repeatCalc == 0 && exist(strcat(Path.output, filesep, 'StimulationPropertiesSUA', ...
                filesep, experiment.name, filesep, 'StimulationProperties_raw','.mat'))
            disp(['Already calculated stimProperties for: ' experiment.name ', expNumber ' num2str(n_animal)])
        else
            cd(Path_stimOld);
            file_ramp=strcat(Path_stimOld,experiment.name,'_Stim',StimRegion,'_ramp.mat');
            file_squa=strcat(Path_stimOld,experiment.name,'_Stim',StimRegion,'_squa.mat');
            
            if exist(file_ramp) && exist(file_squa)
                StimulationProperties_raw=[];
                file_ramp=load( file_ramp);
                file_squa=load( file_squa);
                StimulationProperties_raw=[file_ramp.StimulationProperties;file_squa.StimulationProperties];
                
                %% Save
                if ~isempty(StimulationProperties_raw)
                    if save_data == 0
                    elseif save_data == 1
                        if ~exist(strcat(Path.output,filesep,'StimulationPropertiesSUA', filesep, experiment.name, filesep))
                            mkdir(strcat(Path.output,filesep,filesep,'StimulationPropertiesSUA', filesep, experiment.name, filesep));
                        end
                        save(strcat(Path.output,  filesep, 'StimulationPropertiesSUA', filesep, experiment.name, filesep, 'StimulationProperties_raw','.mat'),'StimulationProperties_raw');
                    end
                else
                    disp('StimProperties was not calculated')
                end
            end
        end
        
    end
end