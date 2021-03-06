function [] = main_function_correctStimData_Stimulation(Experiment,exp_num,Path)

%written by Sabine Gretenkord, last modified 18.07.2017

for iExperiment=1:length(exp_num)
    
    f=exp_num(iExperiment); % f as in Experiment list
    fileName=Experiment(iExperiment).name;
    fileName_forFigTit=strrep(Experiment(iExperiment).name,'_',':');
    
    %get stim properties optogenetics
    optoStimPropFile_folder=strcat(Path.temp,filesep,'Opto_stimProperties');
    optoStimPropFile_file=strcat(optoStimPropFile_folder,filesep,fileName,'.mat');
    load(optoStimPropFile_file)
    
    stimDur=cell2mat(StimulationProperties(:,3));
    stimStart=cell2mat(StimulationProperties(:,11)); %in seconds, from 0
    stimEnd=cell2mat(StimulationProperties(:,12)); %in seconds, from 0
    
    expectedStimTypes=[5,12,900];
    expextedCounteachStim1=[90,30,1];
    expextedCounteachStim2=[90,30,0];
    stimType=expectedStimTypes;
    %stimType=unique(stimDur)';
    %does not work setdiff(expectedStimTypes,stimType)
    
    n_stimType=length(expectedStimTypes);
    prepost_time=stimType;
    
    
    for iStimType=1:n_stimType;
        Logi_StimbyType{iStimType}=stimDur==stimType(iStimType);
        n_StimbyType(iStimType)=sum( Logi_StimbyType{iStimType});
    end
    
    Logi_StimbyType_OBStimulation= Logi_StimbyType;
    
    if   sum(n_StimbyType == expextedCounteachStim1)==3
        
        disp(f)
        disp('correct, 3 Types of stimulations')
        disp(n_StimbyType) 
        
    elseif  sum(n_StimbyType == expextedCounteachStim2)==3
        
        disp(f)
        disp('correct, 2 Types of stimulations')
        disp(n_StimbyType) 
    
    else

        disp(f)
        disp('check stimulation times!')
        disp(n_StimbyType) 
        
    end
    
    
end %iExperiment


end





