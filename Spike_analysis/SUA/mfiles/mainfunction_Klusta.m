function mainfunction_Klusta(experiments,channels,probetype,PRMfolder,DATfolder,BATfolder,Call_path_Anaconda,flag_removeErrChan,flag_T)
% This is a master script that runs the spike sorting algorithm klusta 
% (https://klusta.readthedocs.io/en/latest/) in a format that is friendly 
% to the neuralynx format and the way we save our recordings
% DATfolder should on Q;
%% load experiments and make animals (instead of experiment) the unit on which to sort

% this part here follows the rationale of having multiple recordings 
% ("experiment" in the jergon of get_experiment_redux) from the same
% animal (i.e. a long opto experiment that was subdivided in
% multiple shorter recordings). This piece of code groups such recordings
% all together so that they can be spike-sorted all at the same time.
% If this is not the case for you, you can rearrange it to run over single
% recordings ("experiments") instead. 
 
try
    animals = extractfield(experiments, 'animal_ID');
catch
    animals = extractfield(experiments, 'name');
end
    
animals = animals(~cellfun('isempty', animals));
animals = unique(cellfun(@num2str, animals, 'un', 0));

%% create PRB, PRM and BAT files

% creates PRB and PRM files and puts them into appropriate folders. Since
% the python shell only runs on files that are on C:/, the PRMfolder has 
% to be on C:/. However, since the output files are quite big, the 
% folder2save (where they end up) should better be on Q:/, unless you are 
% just trying out stuff and have space on your local disk. In this case, 
% having the files on your disk will significantly speed up the whole 
% process.

% The PRB file is a file that describes the geometry of the probe

% The PRM file is a file that specifies your preferences for the spike
% sorting process and where to find all the relevant files

% The DAT file is the format in which the recording will be read by the 
% sorting algorithm

% select the channels that you want to spike sort (as they are named in the
% neuralynx format)

% select the probe that you used for your experiments. probe currently
% available are '4shank' and '16_50m'

for idx_animal = 1 : length(experiments)
    disp(strcat('writing animal n°_', num2str(idx_animal)))
    experiment = experiments(idx_animal);
    
    if strcmp(channels, 'HP')
        channels=[experiment.HPreversal-2:experiment.HPreversal+2];
    end 
    
    if strcmp(channels, 'PFC2')
        channels=17:20;
    end
    
    if strcmp(channels, 'PFC5')
        channels=29:32;
    end 
        
    
    try
        animal = experiment.animal_ID;
    catch
        animal = experiment.name;
    end
           
    DATfolder_animal = strcat(DATfolder, animal,'_Chan',num2str(channels(1)),'_',num2str(channels(end))); % folder in which your .DAT and .PRB files will be saved 
                                                                                      % (better be on Q, unless you have a lot of space                                                                                       % on your disk)
    PRB2folder(probetype, DATfolder_animal) % copies PRB file in DATfolder
    

    % create DAT from neuralynx files
    nlx2DAT(animal, experiments,channels, DATfolder_animal,flag_removeErrChan,flag_T)
    
    PRMfolder_animal = strcat(PRMfolder, animal,'_Chan',num2str(channels(1)),'_',num2str(channels(end)));
    PRM2folder(animal, DATfolder_animal, PRMfolder_animal);
end

%% create BAT file 
% A BAT file is a script in lynux. Running klusta with a BAT file permits
% you to run it in batch mode, without having to go through a python
% shell.
createBAT(animals, BATfolder, PRMfolder,channels,Call_path_Anaconda); % CHANGE ONE PATH INSIDE!!!



