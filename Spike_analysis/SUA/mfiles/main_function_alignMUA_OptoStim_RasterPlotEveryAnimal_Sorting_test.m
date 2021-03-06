function main_function_alignMUA_OptoStim_RasterPlotEveryAnimal_Sorting_combineChannles(Experiment,Path,expectedStimTypes,Flag_Response)
%written by Xiaxia, 2019.11.23
Pre_Post_StimPeriod=3;
for iExperiment=1:length(Experiment)
    iExperiment
    fileName=Experiment(iExperiment).animal_ID;
    FileName=Experiment(iExperiment).name;
    StimRegion=Experiment(iExperiment).StimRegion;
    
    chan=Experiment(iExperiment).HPreversal;
    
    % to get the real time 
    File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC1.ncs');
    [T, ~, ~] = load_nlx(File,[]);
    
    t_trueStart=T.TimeStamps_true(1);
        
   if strcmp(Flag_Response,'HP')
       ChanS=[ chan-2, chan-1, chan, chan+1, chan+2];
       fileR=strcat(Path.output,filesep,'SpikeSorting',filesep,'Chan1_16\');
   elseif strcmp(Flag_Response,'PFC2')
       ChanS=17:24;
       fileR=strcat(Path.output,filesep,'SpikeSorting',filesep,'Chan1_16\');
   elseif strcmp(Flag_Response,'PFC5')
       ChanS=25:32;
       fileR=strcat(Path.output,filesep,'SpikeSorting',filesep,'Chan1_16\');
   end
   
    if ~exist(  strcat(strcat(Path.output,filesep,'AlignMUA_RasterPlotEveryAnimal_SpikeSorting',filesep))  )
        mkdir( strcat(strcat(Path.output,filesep,'AlignMUA_RasterPlotEveryAnimal_SpikeSorting',filesep)) )
    end
    
    if exist(  strcat(fileR,fileName,'.mat') ) 
       load( strcat(fileR,fileName,'.mat') );
       if ~isempty(SUAinfo{1,1})
           TimePoint=[];Chan=[];
           for ichan=1:length(SUAinfo{1,1})
               if SUAinfo{1,1}(ichan).channel>ChanS(1)-1 && SUAinfo{1,1}(ichan).channel<ChanS(end)+1
                   TimePoint=[TimePoint;SUAinfo{1, 1}(ichan).Timestamps] ;
               end
           end
           
           TimePoint=unique(TimePoint); % SUA, point
           TimePoint=1000*TimePoint/32000+t_trueStart/1000;
           
           % load stimulation
           optoStimPropFile_folder=strcat(Path.temp,filesep,'Opto_stimProperties');
           if ismember(StimRegion,{'HP','PFC5','PFC2','LEC'})
               optoStimPropFile_file=strcat(optoStimPropFile_folder,filesep,FileName,'.mat');
               load(optoStimPropFile_file)
           elseif strcmp(StimRegion,'PFC5+PFC2')
               optoStimPropFile_file=strcat(optoStimPropFile_folder,filesep,FileName,'_Stim',Expect_Region,'.mat');
               load(optoStimPropFile_file)
           end
           
           stimTypes=cellstr(StimulationProperties(:,8));
           stimType=unique(stimTypes)';
           
           n_stimType=length(expectedStimTypes);
           for iStimType=1:n_stimType;
               StimType=expectedStimTypes{iStimType};
               if sum (strcmp(StimType, stimType))
                   
                   Logi_StimbyType=strcmp(stimTypes,expectedStimTypes(iStimType));
                   n_StimbyType=sum( Logi_StimbyType);
                   index_StimbyType=find (Logi_StimbyType~=0);
                   StimulationProperties_StimbyType=StimulationProperties(index_StimbyType,:);
                   num=0; TempSpike=[];
                   for num_stim=1:n_StimbyType
                       StimulationProperties_StimbyType_each=StimulationProperties_StimbyType(num_stim,:);
                       stimPeriod=StimulationProperties_StimbyType_each{3};
                       if strcmp(expectedStimTypes{iStimType}(1:4),'squa')&& stimPeriod==2
                           E=cell2mat( StimulationProperties_StimbyType_each(1,2) )+(Pre_Post_StimPeriod+1)*10^3; % ms
                       end
                       E=cell2mat( StimulationProperties_StimbyType_each(1,2) )+Pre_Post_StimPeriod*10^3; % ms
                       S=cell2mat( StimulationProperties_StimbyType_each(1,1) )-Pre_Post_StimPeriod*10^3; % ms
                       index_s=min(find(TimePoint>S));
                       index_e=max(find(TimePoint<E));
                       if index_e>index_s
                           num=num+1;
                           index_spike=((TimePoint(index_s:index_e)-S)/1000);
                           TempSpike{num,1}=index_spike';
                       end
                   end
                   
                  
                   if ~isempty( TempSpike )
                       MarkerFormat.MarkerSize = 4;
                       plotSpikeRaster(TempSpike,'PlotType','scatter','XLimForCell',[0 9],'MarkerFormat',MarkerFormat);
                       title (strcat (StimType));
                       xlim([2 7])
                       colormap('jet')
                       pathout=strcat(Path.output,filesep,'AlignMUA_RasterPlotEveryAnimal_SpikeSorting',filesep);
                       saveas(figure(1),strcat(pathout,FileName,'_',StimType,'_Stim_',StimRegion,'_Res_',Flag_Response,'.eps'));
                       saveas(figure(1),strcat(pathout,FileName,'_',StimType,'_Stim_',StimRegion,'_Res_',Flag_Response,'.jpg'));
                       close(figure (1));
                   end
                   
               end
               
           end
           
       end
       
    end
    
end
                   
                   
                       
                 







