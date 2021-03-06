function [ output_args ] = main_function_loadLFP(Experiment,Path)
%written by Sabine Gretenkord - last modified 24.10.2016

%Loads Neuralynx (.ncs) data: two selected channels (PL and HP),
%Downsamples after appropriate filtering (Nyquist frequency)
%Filters the signal so that the output can be considered LFP data (<300 Hz)
%Saves data in the 'Temp' folder

%required functions: load_nlx,filter_downsample,ZeroPhaseFilterZeroPadding
%found here: Q:\directionality meeting\Scripts\LoadFilterDownsample

for iExperiment=1:length(Experiment)
    
    Experiment(iExperiment).chan_num{1}= Experiment(iExperiment).PLchoice;
    Experiment(iExperiment).chan_num{2}= Experiment(iExperiment).HPchoice;
    
    for iChan=1:2
        if Experiment(iExperiment).chan_num{iChan}~=0
            
            %load data (amplitude adjustment automatically)
            thisChannel=Experiment(iExperiment).chan_num{iChan};
            ExtractModeArray=[];
            File= strcat(Experiment(iExperiment).path,Experiment(iExperiment).name,'\CSC',num2str(thisChannel),'.ncs');
            [TimeStamps, Samples, fs] = load_nlx(File,ExtractModeArray);
            
            %downsample
            DownsampleFactor=32; %Downsample to sampling rate fs: 1000 Hz, after filtering high cut: 500 Hz
            [TimeStamps, Samples, fs] = filter_downsample(TimeStamps, Samples, fs, DownsampleFactor);
            
            %filter for LFP band (<300 Hz)
            high_cut=300; %LFP high cut filter
            %Samples=ZeroPhaseFilterZeroPadding(Samples,fs,[0 high_cut]);
            Samples=ZeroPhaseFilterZeroPadding(Samples,fs,[30 48]);
            %save data in temp folder
            OutputFolderName=strcat('nlx_load_LFP');
            if ~exist(  strcat(Path.temp,filesep,OutputFolderName,filesep,Experiment(iExperiment).name)  )
                mkdir(  strcat(Path.temp,filesep,OutputFolderName,filesep,Experiment(iExperiment).name)  )
            end
            save(strcat(Path.temp,filesep,OutputFolderName,filesep,Experiment(iExperiment).name,...
                filesep,'CSC',num2str(thisChannel)),...
                'TimeStamps','Samples','fs','high_cut','DownsampleFactor')
        end
    end

end

end

