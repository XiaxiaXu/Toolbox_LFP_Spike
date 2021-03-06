
%% load experiment data and organize by animals
clear all
load('E:\Projects\Project-LEC\Experiment-InjectionLEC\Experiment list\Experiment_injectionLEC_StimHP_Con.mat')

%% process & save timestamps and waveform features
StimRegion='HP';ResponseRegion='HP';Chans=[1 16];
StimRegion='PFC5';ResponseRegion='PFC';Chans=[17 32];
StimRegion='PFC2';ResponseRegion='PFC';Chans=[17 32];
StimRegion='LEC';ResponseRegion='HP';Chans=[1 16];
StimRegion='LEC';ResponseRegion='PFC';Chans=[17 32];

Path.SUA=strcat('F:\Klusta\Stim',StimRegion);
PRMpath = strcat(Path.SUA,'\PRMfolder\'); % main folder in which you store your PRM files
DATpath=strcat(Path.SUA,'\DATfolder\');

save_data = 1;         %save data or not
directory2save = strcat(Path.SUA,'\SortedSUA'); % where you want to save the output

for iExperiment=1:length(Experiment)
    disp(strcat('working on animal n°:', num2str(iExperiment)))
    animal =Experiment(iExperiment).animal_ID;
    DATfolder =strcat(DATpath,animal,'_Chan',num2str(Chans(1)),'_',num2str(Chans(2)));% folder in which your .DAT and .PRB files are saved
    PRMfolder =strcat(PRMpath,animal,'_Chan',num2str(Chans(1)),'_',num2str(Chans(2)));% folder in which your .DAT and .PRB files are saved
    if exist(DATfolder) && exist(PRMfolder)
        loadKlusta(animal, PRMfolder, DATfolder, save_data, directory2save);
    end
end

% seperate the recording files

for iExperiment=1%:length(Experiment)
    disp(strcat('working on animal n°:', num2str(iExperiment)))
    animal =Experiment(iExperiment).animal_ID;
    load(strcat(Path.SUA,'\SortedSUA\', filesep, animal,'.mat'))
    [r,c]=size(SUAinfo(1,:));
    for num=1:1
        SUA_struct=SUAinfo{1,num};
        fileName=Experiment(iExperiment).name;
%         fileName=SUA_struct.file;
        disp(strcat('working on file:', fileName))
        for unit=1: length( SUA_struct)
            unit
            SUA = strcat('SUA', num2str(unit));
            SUAstruct.(SUA).spiketimes = SUA_struct(unit).Timestamps'/32000;
            SUAstruct.(SUA).recording_site = 2;
            [~,locAmp]=min(SUA_struct(unit).Waveform');
            SUAstruct.(SUA).waveform= SUA_struct(unit).Waveform(locAmp-21:locAmp+30)';
            SUAstruct.(SUA).CSC = SUA_struct(unit).channel+16;
            SUAstruct.(SUA).num_spikes = length(SUA_struct(unit).Timestamps);
            SUAstruct.(SUA).ISI=diff(SUA_struct(unit).Timestamps)/32000;
            SUAstruct.(SUA).cluster=1;
            SUAstruct.(SUA).WaveFeatures=SUA_struct(unit).MeanWaveFeatures;
            
            if ~exist(strcat(directory2save, filesep, fileName))
                mkdir(strcat(directory2save, filesep, fileName));
            end
        end
        if save_data == 1
            save(strcat(directory2save, filesep, fileName,filesep,'SUAstruct'), 'SUAstruct')
        end
        clear SUA_struct fileName SUAstruct SUA unit
    end
end

%% getStimProperties for SUA i.e. 
repeatCalc=1;
SUAgetStimulationProperties(Experiment, save_data, repeatCalc,Path)
%SUAcorrectStimulationProperties(Experiment, save_data, repeatCalc,Path)
for i=1:length(Experiment)
    experiment=Experiment(i);
    experiment.animal_ID
    SUAdata = getStimulationSUAspikesKlusta(experiment, save_data, repeatCalc,Path);
end

