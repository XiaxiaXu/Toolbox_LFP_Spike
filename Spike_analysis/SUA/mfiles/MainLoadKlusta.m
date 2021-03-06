%% load experiment data and organize by animals

clear all
klusta = 2;
experiments = get_experiment_redux(klusta);
% experiments = experiments(903 : end);
experiments = experiments(~cellfun(@isempty,{experiments.animal_ID}));
% this part here follows the rationale of having multiple recordings 
% ("experiment" in the jergon of get_experiment_redux) from the same
% animal (i.e. a long opto experiment that was subdivided in
% multiple shorter recordings). This piece of code groups such recordings
% all together so that they can be spike-sorted all at the same time.
% If this is not the case for you, you can rearrange it to run over single
% recordings ("experiments") instead. 


animals = extractfield(experiments, 'animal_ID');
animals = animals(~cellfun('isempty', animals));
[animals, unique_animals_idx] = unique(cellfun(@num2str, animals, 'un', 0));

%% process & save timestamps and waveform features

PRMfolder = 'E:\klusta\PRM files\'; % main folder in which you store your PRM files
                                    % (which then are in an animal-specific subfolder)
save_data = 1;
directory2save = 'Q:\Personal\Mattia\PFCmicro\results\SUA Klusta'; % where you want to save the output
for idx_animal = 1 : length(animals)
    disp(strcat('working on animal n°_', num2str(idx_animal)))
    animal = animals{idx_animal};
    DATfolder = strcat('Q:/Personal/Mattia/PFCmicro/results/klusta/', animal); % folder in which your .DAT and .PRB files are saved
    SUAinfo = loadKlusta(animal, PRMfolder, DATfolder, save_data, directory2save);
end

