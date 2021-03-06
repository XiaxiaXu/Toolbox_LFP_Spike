function main_LoadKlusta_SUAinfo_SUAstruct(Experiment,Path,PRMpath,DATpath,Chans,save_data,repeatCalc,directory2save)

if Chans(2)==16
    recording_site=1;
    Chancorrect=0;
elseif Chans(2)==32
    recording_site=2;
    Chancorrect=16;
elseif Chans(2)==48
    recording_site=3;
    Chancorrect=32;
end

% for iExperiment=1:length(Experiment)
%     animal =Experiment(iExperiment).animal_ID;
%     DATfolder =strcat(DATpath,animal,'_Chan',num2str(Chans(1)),'_',num2str(Chans(2)));% folder in which your .DAT and .PRB files are saved
%     PRMfolder =strcat(PRMpath,animal,'_Chan',num2str(Chans(1)),'_',num2str(Chans(2)));% folder in which your .DAT and .PRB files are saved
%     if exist(DATfolder) && exist(PRMfolder)
%         disp(strcat('working on animal n°:', num2str(iExperiment)))
%         loadKlusta(animal, PRMfolder, DATfolder, save_data, directory2save);
%     end
% end

% seperate the recording files
for iExperiment=1:length(Experiment)
    disp(strcat('working on animal n°:', num2str(iExperiment)))
    animal =Experiment(iExperiment).animal_ID;
    if exist(strcat(Path.SUA,'\SortedSUA\', filesep, animal,'.mat'))
        load(strcat(Path.SUA,'\SortedSUA\', filesep, animal,'.mat'))
        [r,c]=size(SUAinfo(1,:));
        for num=1:c
            SUA_struct=SUAinfo{1,num};
            fileName=Experiment(iExperiment).name;
            %         fileName=SUA_struct.file;
            disp(strcat('working on file:', fileName))
            for unit=1: length( SUA_struct)
                SUA = strcat('SUA', num2str(unit));
                SUAstruct.(SUA).spiketimes = SUA_struct(unit).Timestamps'/32000;
                SUAstruct.(SUA).recording_site = recording_site;
                [~,locAmp]=min(SUA_struct(unit).Waveform');
                
                try
                    SUAstruct.(SUA).waveform=SUA_struct(unit).Waveform(locAmp-21:locAmp+30)';
                    SUAstruct.(SUA).cluster=SUA_struct(unit).ClusterType;
                catch
                    SUAstruct.(SUA).waveform= SUA_struct(unit).Waveform';
                    SUAstruct.(SUA).cluster=0;
                end
                
                SUAstruct.(SUA).CSC = SUA_struct(unit).channel+Chancorrect;
                SUAstruct.(SUA).num_spikes = length(SUA_struct(unit).Timestamps);
                SUAstruct.(SUA).ISI=diff(SUA_struct(unit).Timestamps)/32000;
                SUAstruct.(SUA).WaveFeatures=SUA_struct(unit).MeanWaveFeatures;
                
                if ~exist(strcat(directory2save, filesep, fileName))
                    mkdir(strcat(directory2save, filesep, fileName));
                end
            end
            
        end
        if save_data == 1
            save(strcat(directory2save, filesep, fileName,filesep,'SUAstruct'), 'SUAstruct')
        end
        clear SUA_struct fileName SUAstruct SUA unit
    end
end



