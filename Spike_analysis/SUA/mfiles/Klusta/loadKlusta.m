function [SUAinfo, features] = loadKlusta(animal, PRMfolder, DATfolder, save_data, directory2save)
% by Mattia 25/01

% this function saves timestamps for every SUA, spike wave features and
% parameters of the cluster (for quality assesment)
% in this version, it loops through consecutive recordings of the same
% animal (this is how my experiments are). if this is not the case for you,
% get read of the most outer loop and adjust the code accordingly.

% input :   - animal: string
%           - PRMfolder: where the PRM file can be found - string
%           - DATfolder: where the DAT file can be found) - string
%           - save_data: 1 yes, 0 no (only loading)
%           - directory2save: where to save the output - string

% output :  - SUAinfo: structure with timestamps and cluster info
%           - features: waveform features for PYR/IN clustering (read
%             getSpikeFeatures documentation for more info)




filenamekwik = strcat(PRMfolder, filesep, animal, '.kwik'); % name kwik file
filenamekwx = strcat(PRMfolder,filesep, animal, '.kwx'); % name kwix file

TimeStamps = double(hdf5read(filenamekwik, '/channel_groups/0/spikes/time_samples')); % timestamps for all detected spikes
RecordingBreaks = find(diff(TimeStamps) < 0); % only useful if you have more than one recording per animal
Clusters = hdf5read(filenamekwik, '/channel_groups/0/spikes/clusters/main'); % assign a cluster to every spike
FeatureMask = hdf5read(filenamekwx, '/channel_groups/0/features_masks'); % load the three features on every channel for every spike

filePattern = fullfile(DATfolder, '*.dat');
DATfiles = dir(filePattern);

for DATfile_idx = 1 : numel(DATfiles)
    
    DATfile = DATfiles(DATfile_idx).name;
    gwfparams.file = strcat(DATfolder, filesep, DATfile);                                 % path to your DAT file
    gwfparams.fileName = DATfile;                                                         % DAT file name
    gwfparams.dataType = 'int16';                                                         % Data type of .dat file (there should be no need ot change this)
    gwfparams.nCh = 16;                                                                   % Number of channels that were clustered
    gwfparams.wfWin = [-120 120];                                                           % Number of samples before and after spiketime to include in waveform
    gwfparams.nWf = 500;                                                                  % Number of waveforms per unit to pull out
    if DATfile_idx == 1
        if numel(RecordingBreaks) > 0
            gwfparams.spikeTimes = TimeStamps(1 : RecordingBreaks(DATfile_idx));              % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = Clusters(1 : RecordingBreaks(DATfile_idx));             % Vector of cluster IDs (Phy nomenclature) same length as .spikeTimes
        else
            gwfparams.spikeTimes = TimeStamps;
            gwfparams.spikeClusters = Clusters;
        end
    elseif DATfile_idx == numel(DATfiles)
        gwfparams.spikeTimes = TimeStamps(RecordingBreaks(DATfile_idx - 1)...
            : end);
        gwfparams.spikeClusters = Clusters(RecordingBreaks(DATfile_idx - 1)...
            : end);
    else
        gwfparams.spikeTimes = TimeStamps(RecordingBreaks(DATfile_idx - 1)...
            : RecordingBreaks(DATfile_idx));
        gwfparams.spikeClusters = Clusters(RecordingBreaks(DATfile_idx - 1)...
            : RecordingBreaks(DATfile_idx));
    end
    
   
    WaveForms = getWaveForms(gwfparams);                                                  % function from CortexLab extracting Waveforms for each cluster
    %WaveForms.waveFormsMean = WaveForms.waveFormsMean - mean(WaveForms.waveFormsMean, 3); % zero-center the mean Waveforms
    SUA = struct;
    countSUA = 0;
    fs = 32000;
    refractory_period = fs * 2 * 10 ^-3; % two ms in timestamps
    
    for idxCluster = 1 : length(WaveForms.unitIDs) % loop over clusters
        
        ClusterType = hdf5read(filenamekwik, strcat( '/channel_groups/0/clusters/main/',  ...
            num2str(WaveForms.unitIDs(idxCluster)), '/cluster_group')); % 0 for noise, 1 for MUA and 2 for SUA
        
        if ClusterType == 2 || ClusterType == 1 % if SUA or mua
            
            countSUA = countSUA + 1;
            [~, SUA(countSUA).channel] = min(min(WaveForms.waveFormsMean(idxCluster, :, :), [], 3)); % channel is the one with biggest negativity
            % in the average mask
            SUA(countSUA).Timestamps = gwfparams.spikeTimes(gwfparams.spikeClusters == ...
                WaveForms.unitIDs(idxCluster));                                                      % timestamps when spike occurs
            SUA(countSUA).ClusterID = WaveForms.unitIDs(idxCluster);                                 % to check stuff in phy, you never know
            SUA(countSUA).Waveform = squeeze(WaveForms.waveFormsMean(idxCluster, ...
                SUA(countSUA).channel, :));                                                          % for further analysis, i.e. distinguishing PYR and IN
            SUA(countSUA).Amplitudes = squeeze(FeatureMask(1, SUA(countSUA).channel * 3 - 1, ...
                Clusters == WaveForms.unitIDs(idxCluster)));                                         % the second feature, on the best channel (see above). for cluster quality assesment
            SUA(countSUA).RPV = nnz(diff(gwfparams.spikeTimes) < refractory_period)/length(gwfparams.spikeTimes);                 % refractory period violations. for cluster quality assesment
            SUA(countSUA).MeanWaveFeatures = getSpikeFeatures(SUA(countSUA).Waveform);               % for SUA clustering
            %features(countSUA, :) = SUA(countSUA).MeanWaveFeatures;
            %             SUA(countSUA).file = strip(strip(strip(strip(DATfile, 'right', 't'), ...                 % stupid way of getting to the name
            %                 'right', 'a'), 'right', 'd'), 'right', '.');                                         % (deleting .DAT at the end).
            SUA(countSUA).file = DATfile;
            SUA(countSUA).ClusterType =ClusterType;
        end
        
        
    end
    SUAinfo{DATfile_idx} = SUA;
end


if save_data == 1
    if ~exist(directory2save)
        mkdir(directory2save)
    end
    save(strcat(directory2save, filesep, animal), 'SUAinfo')
end