function main_function_alignMUA_OptoStim_RasterPlotEveryAnimal_Sorting(Experiment,Path,expectedStimTypes,Flag_Response)
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
    
    t_trueStart=T.TimeStamps_true(1)/10^3; % ms
    
    if strcmp(Flag_Response,'HP')
%         ChanS=[ chan-2, chan-1, chan, chan+1, chan+2];
        ChanS=1:16;
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
            TimePoints=[];num_cell=0;
            for ichan=1:length(SUAinfo{1,1})
                if SUAinfo{1,1}(ichan).channel>ChanS(1)-1 && SUAinfo{1,1}(ichan).channel<ChanS(end)+1
                    num_cell=num_cell+1;
                    TimePoints{num_cell}=1000*SUAinfo{1, 1}(ichan).Timestamps/32000-t_trueStart; % ms
                end
            end
            
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
                    
                    m=ceil(sqrt(length(TimePoints)));
                    n=m;
                    
                    for num_cell=1:length(TimePoints)
                        TimePoint=TimePoints{num_cell};
                        
                        num=0; TempSpike=[];Spikes=[];
                        for num_stim=1:n_StimbyType
                            StimulationProperties_StimbyType_each=StimulationProperties_StimbyType(num_stim,:);
                            stimPeriod=StimulationProperties_StimbyType_each{3};
                            if strcmp(expectedStimTypes{iStimType}(1:4),'squa')&& stimPeriod==2
                                E=cell2mat( StimulationProperties_StimbyType_each(1,2) )+(Pre_Post_StimPeriod+1)*10^3; % ms
                            end
                            E=cell2mat( StimulationProperties_StimbyType_each(1,2) )+Pre_Post_StimPeriod*10^3; % ms
                            S=cell2mat( StimulationProperties_StimbyType_each(1,1) )-Pre_Post_StimPeriod*10^3; % ms
                            
%                             File= strcat(Experiment(iExperiment).path,filesep,Experiment(iExperiment).name,'\CSC',num2str(SUAinfo{1,1}(num_cell).channel),'.ncs');
%                             [~, recordingRaw, samplingrate_MUA] = load_nlx(File,[S E]);
%                             recordingMUA = ZeroPhaseFilter(recordingRaw,samplingrate_MUA,[500 5000]);
%                             thr = std( recordingMUA)*threshold;
%                             [peakLoc, ~] = peakfinderOpto( recordingMUA,thr/2 ,-thr,-1,false);
%                             
                            index_s=min(find(TimePoint>=S));
                            index_e=max(find(TimePoint<=E));
                            
                            if index_e>index_s
                                num=num+1;
                                index_spike=((TimePoint(index_s:index_e)-S)/1000);
                                TempSpike{num,1}=index_spike';
%                                 spikes=zeros(1,10*32000);
%                                 ispike=round(index_spike*32000);
%                                 spikes(ispike)=1;
%                                 Spikes=[Spikes;spikes];
                            end
                            
                        end
                        
                        if ~isempty( TempSpike )
                            subplot(m,n,num_cell);
                            MarkerFormat.MarkerSize = 4;
                            plotSpikeRaster(TempSpike,'PlotType','scatter','XLimForCell',[0 9],'MarkerFormat',MarkerFormat);
                            title (strcat (StimType));
                            xlim([2 7])
                            colormap('jet')
                        end
                    end
                    pathout=strcat(Path.output,filesep,'AlignMUA_RasterPlotEveryAnimal_SpikeSorting',filesep);
                    saveas(figure(1),strcat(pathout,FileName,'_',StimType,'_Stim',StimRegion,'_Res',Flag_Response,'.eps'));
                    saveas(figure(1),strcat(pathout,FileName,'_',StimType,'_Stim',StimRegion,'_Res',Flag_Response,'.jpg'));
                    close(figure (1));
                    
                end
                
            end
            
        end
        
    end
    
end











