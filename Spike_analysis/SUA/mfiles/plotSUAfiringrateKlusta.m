function spikes_tot = plotSUAfiringrateKlusta(Experiment, expType, StimRegion,stimulusType, RespArea, age, cluster_idx, repeatCalc)
%% By Mattia
Path = get_path;
Gwindow = gausswin(1001, 10); % gaussian window of 1000ms with stdev of 100ms
Gwindow = Gwindow / sum(Gwindow); % normalize the gaussian kernel
save_data = 1;
spikes_tot = [];
firing_tot = [];
OMI = [];
OMI2 = [];
pvalue = [];
pvalue2 = [];
ISIpre = [];
ISIduring = [];
ISIpost = [];
corr_coeffs_stim = [];
corr_coeffs_pre_stim = [];
corr_coeffs_post_stim = [];
pop_coupling = [];
channels = [];
pre_stim = 1 : 3000; % in ms, ramp format
stim = 3001 : 6000; % in ms, ramp format
post_stim = 6001 : 9000; % in ms, ramp format
stim_layer = [];

for n_animal =1: length(Experiment)
    experiment = Experiment(n_animal);
    if nnz( ~ isnan(experiment.StimRegion)) && strcmp(expType,experiment.Exp_type)
        if ~isempty(experiment.animal_ID)  && ismember(experiment.StimRegion, StimRegion) && ismember(experiment.AgeP, age)
            if strcmp(RespArea, 'HP')
                CSC = experiment.HPreversal - 2 : experiment.HPreversal + 2;
            elseif strcmp(RespArea, 'PFC2')
                CSC = 17 : 24;
            elseif strcmp(RespArea, 'PFC5')
                CSC = 25 : 32;
            elseif strcmp(RespArea, 'PFC')
                CSC = 17 : 32;
            end
            SUAdata = getStimulationSUAspikesKlusta(experiment, save_data, repeatCalc,Path);
            CSCs = SUAdata.channels + 16;
            if strcmp(stimulusType,'ramp') % square ramp
            spikes_animals = SUAdata.ramp_spike_matrix;  
            elseif strcmp(stimulusType,'square') % square ramp
            spikes_animals = SUAdata.pulse_spike_matrix;
            end
            clusters = SUAdata.clusters;
            if cluster_idx > 0 %% collect the right channel, remove some cluster if need; cluster_idx=0 means use all
                spikes_animals(:, clusters ~= cluster_idx | ~ismember(CSCs, CSC), :) = [];
                CSCs(clusters ~= cluster_idx) = [];
            else
                spikes_animals(:, ~ismember(CSCs, CSC), :) = [];
            end
            spikes_animal=[];
            if numel(spikes_animals) > 0 && size(spikes_animals, 2) > 1
                    iClusters=0;
                    for iC=1:size(spikes_animals, 2)
                        if sum( sum (spikes_animals(:,iC,1:3000)))>10
                            iClusters=iClusters+1;
                         spikes_animal(:,iClusters,:) = spikes_animals(:,iC,:);
                        end
                    end
            end
             if numel(spikes_animal) > 0 
            % first concatenate the spike tensor into a matrix
            spikes_convolved = reshape(permute(spikes_animal, [2 3 1]), size(spikes_animal, 2), []);
            % convolve it with a gaussian window for better corr estimation
            
            for unit = 1 : size(spikes_convolved, 1)
                spikes_convolved(unit, :) = conv(spikes_convolved(unit, :), Gwindow, 'same');
            end
            
            % reshape it back so that you have separated trials of stimulation
            spikes_convolved= permute(reshape(spikes_convolved, size(spikes_animal, 2),size(spikes_animal, 3), []), [3 1 2]);
            goodcluster=0;
           for  cluster = 1 : size(spikes_convolved, 1)
               if sum(spikes_convolved(cluster,:))>0
                goodcluster=goodcluster+1;
                spikes_convolved(goodcluster,:)=spikes_convolved(cluster,:);
               end
           end
            % initialize variables as NaN
            corr_coeffs_animal_stim = NaN(size(spikes_convolved, 1), size(spikes_convolved, 2), size(spikes_convolved, 2));
            corr_coeffs_animal_pre_stim = corr_coeffs_animal_stim ;
            corr_coeffs_animal_post_stim = corr_coeffs_animal_stim;
            % loop over trials and collect corr coeffs
            % fisher-transform (or not) the pearson coefficients
            for stims = 1 : size(spikes_convolved, 1)
                % THIS IS WITH FISHER CORRECTION
%                 corr_coeffs_animal_stim(ramp, :, :) = atanh(corrcoef(squeeze(spikes_convolved(ramp, :, stim))'));
%                 corr_coeffs_animal_pre_stim(ramp, :, :) = atanh(corrcoef(squeeze(spikes_convolved(ramp, :, pre_stim))'));
%                 corr_coeffs_animal_post_stim(ramp, :, :) = atanh(corrcoef(squeeze(spikes_convolved(ramp, :, post_stim))'));
                % THIS IS WITHOUT
                corr_coeffs_animal_stim(stims, :, :) = corrcoef(squeeze(spikes_convolved(stims, :, stim))');
                corr_coeffs_animal_pre_stim(stims, :, :) = corrcoef(squeeze(spikes_convolved(stims, :, pre_stim))');
                corr_coeffs_animal_post_stim(stims, :, :) = corrcoef(squeeze(spikes_convolved(stims, :, post_stim))');
            end
            % average over stims
            corr_coeffs_animal_stim = nanmean(corr_coeffs_animal_stim);
            corr_coeffs_animal_pre_stim = nanmean(corr_coeffs_animal_pre_stim);
            corr_coeffs_animal_post_stim = nanmean(corr_coeffs_animal_post_stim);
            % concatenate in arrays with corr_coefs for all animals
            corr_coeffs_stim = cat(1, corr_coeffs_stim, corr_coeffs_animal_stim(:));
            corr_coeffs_pre_stim = cat(1, corr_coeffs_pre_stim, corr_coeffs_animal_pre_stim(:));
            corr_coeffs_post_stim = cat(1, corr_coeffs_post_stim, corr_coeffs_animal_post_stim(:));
            
          
                if size(spikes_animal, 2) > 1
                    spikes_units = squeeze(mean(spikes_animal));
                else
                    spikes_units = squeeze(mean(spikes_animal))';
                end
                firing_units(:, 1) = log10(mean(spikes_units(:, pre_stim), 2));
                firing_units(:, 2) = log10(mean(spikes_units(:, stim), 2));
                firing_units(:, 3) = log10(mean(spikes_units(:, post_stim), 2));
                spikes_tot = cat(1, spikes_tot, spikes_units);
                firing_tot = cat(1, firing_tot, firing_units);
                pre = squeeze(sum(spikes_animal(:, :, pre_stim), 3));
                during = squeeze(sum(spikes_animal(:, :, stim), 3));
                post = squeeze(sum(spikes_animal(:, :, post_stim), 3));
                OMI_animal = nanmean((during - pre) ./ (during + pre)); % compute modulation index
                OMI_animal2 = nanmean((post - during ) ./ (during + post)); % compute modulation index
                pvalue_animal = zeros(1, size(pre, 2)); % preallocate
                 pvalue_animal2 = zeros(1, size(pre, 2)); % preallocate
                ISIbins = linspace(0, 200, 50);
                animalISIpre = zeros(size(pre, 2), numel(ISIbins) - 1);
                animalISIduring = animalISIpre; 
                animalISIpost = animalISIpre;
                for unit = 1 : size(pre, 2)% pvalue for each unit
                    pvalue_animal(unit) = signrank(pre(:, unit), during(:, unit)); % compute pvalue of "modulation index"
                     pvalue_animal2(unit) = signrank(post(:, unit), during(:, unit)); % compute pvalue of "modulation index"
                    spiketimes_pre = find(spikes_animal(:, unit, pre_stim));
                    spiketimes_during = find(spikes_animal(:, unit, stim));
                    spiketimes_post = find(spikes_animal(:, unit, post_stim));
                    animalISIpre(unit, :) = histcounts(diff(spiketimes_pre), ISIbins);
                    animalISIduring(unit, :) = histcounts(diff(spiketimes_during), ISIbins);
                    animalISIpost(unit, :) = histcounts(diff(spiketimes_post), ISIbins);
                    
                end
                if nnz(strcmp(expType, 'niet')) && nanmin(pvalue_animal) < 0.01
                    disp(strcat('animal_', experiment.animal_ID, '_has_', ...
                        num2str(nnz(pvalue_animal < 0.01)), '_units with pvalue < 0.01'))
                end
                OMI = horzcat(OMI, OMI_animal); % concatenate
                OMI2 = horzcat(OMI2, OMI_animal2);
                pvalue = horzcat(pvalue, pvalue_animal); % concatenate
                pvalue2 = horzcat(pvalue2, pvalue_animal2); % concatenate
                ISIpre = vertcat(ISIpre, animalISIpre);
                ISIduring = vertcat(ISIduring, animalISIduring);
                ISIpost = vertcat(ISIpost, animalISIpost);
                channels = horzcat(channels, CSCs);
                clear firing_units
                % load pop coupling stuff
%                 load(strcat('Q:\Personal\Mattia\PFCmicro\results\population coupling\', experiment.animal_ID))
%                 if length(OMI_animal) ~= length(pop_coupling_stuff.pop_coupling)
%                     pop_coupling = cat(1, pop_coupling, NaN(length(OMI_animal), 1));
%                 else
%                     pop_coupling = cat(1, pop_coupling, pop_coupling_stuff.pop_coupling);
%                 end
            end
            if strcmp(experiment.StimRegion, 'PFC2')
                stim_layer = cat(1, stim_layer, ones(length(OMI_animal), 1));
            elseif strcmp(experiment.StimRegion, 'PFC5')
                stim_layer = cat(1, stim_layer, ones(length(OMI_animal), 1) * 2);
            end
            clearvars spikes_animal
        end
    end
end

figure 
subplot(2,3,1)
cdfplot(corr_coeffs_stim)
hold on
cdfplot(corr_coeffs_pre_stim)
hold on
cdfplot(corr_coeffs_post_stim)
ylabel('Cumulative proportion'); xlabel('Pearson R'); title('Corr Coeffs pre, stim and post')
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')

% volcano plot; size is equal to number of spikes and colored upon
% significant or not
RdBu = cbrewer('div', 'RdBu', 100);
num_spikes = sum(spikes_tot, 2); % for size only
num_spikes(num_spikes == 0) = 0.1; % scatter can't handle 0 as an input for size
max_spikes = max(num_spikes);
scaling_factor = 200 / max_spikes;

hold on
subplot(2,3,2)
scatter(OMI, - log10(pvalue), num_spikes * scaling_factor, pvalue < 0.01, 'filled');
colormap(flip(RdBu)); caxis([-0.1 1.1])
hold on
plot([0 0], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for no modulation
xlim([-1.1 1.1]); ylim([-0.1 ceil(max(- log10(pvalue)))])
ylabel('- Log10 (pvalue)'); xlabel('Modulation Index stim/pre'); alpha(0.5)
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')
subplot(2,3,3)
scatter(OMI2, - log10(pvalue), num_spikes * scaling_factor, pvalue < 0.01, 'filled');
colormap(flip(RdBu)); caxis([-0.1 1.1])
hold on
plot([0 0], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for no modulation
xlim([-1.1 1.1]); ylim([-0.1 ceil(max(- log10(pvalue)))])
ylabel('- Log10 (pvalue)'); xlabel('Modulation Index post/stim'); alpha(0.5)
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')

% modulation = [nnz(OMI < 0 & pvalue < 0.01) nnz(pvalue > 0.01) nnz(OMI > 0 & pvalue < 0.01)];
modulation = [nnz(OMI < 0 & pvalue < 0.01) nnz(pvalue > 0.01) nnz(OMI > 0 & pvalue < 0.01)];%% mattia use p<0.01
modulation = modulation ./ length(OMI);
subplot(2,3,4)
bar([modulation; nan(1, length(modulation))], 'stacked')
set(gca,'xtick',1); xlim([0 2]);
ylabel('Proportion'); xlabel(''); xticks(''); xticklabels('')
title('Proportion of (un)modulated units stim/pre')
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')

modulation2 = [nnz(OMI2 < 0 & pvalue2 < 0.01) nnz(pvalue2 > 0.01) nnz(OMI2 > 0 & pvalue2 < 0.01)];%% mattia use p<0.01
modulation2 = modulation2 ./ length(OMI2);
subplot(2,3,5)
bar([modulation2; nan(1, length(modulation2))], 'stacked')
set(gca,'xtick',1); xlim([0 2]);
ylabel('Proportion'); xlabel(''); xticks(''); xticklabels('')
title('Proportion of (un)modulated units post/stim')
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')
%% plot normalized ISI __not check yet#
ISIpre = ISIpre(OMI > 0 & pvalue < 0.01, :);
ISIduring = ISIduring(OMI > 0 & pvalue < 0.01, :);

cmap = cbrewer('div', 'Spectral', 100); % fancy colormap
figure;
subplot(1,2,1)
boundedline(linspace(0, 200, size(ISIpre, 2)), mean(ISIpre) ./ sum(mean(ISIpre)), std(ISIpre) ./ ...
    sqrt(size(ISIpre, 1)) ./ sum(mean(ISIpre)), 'cmap', cmap(15, :))
boundedline(linspace(0, 200, size(ISIpre, 2)), mean(ISIduring) ./ sum(mean(ISIduring)), std(ISIduring) ./ ...
    sqrt(size(ISIpre, 1)) ./ sum(mean(ISIduring)), 'cmap', cmap(75, :))
% boundedline(linspace(0, 200, 19), mean(ISIpost) ./ sum(mean(ISIpost)), std(ISIpost) ./ ...
%     sqrt(size(ISIpre, 1)) ./ sum(mean(ISIpost)), 'cmap', cmap(100, :))
ylabel('ISI occurrences'); xlabel('time (ms)')
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')

% spikes_tot_increase=spikes_tot(OMI > 0 & pvalue < 0.01, :);
% spikes_tot_nor=spikes_tot_increase./mean(spikes_tot_increase(:,1:3000),2);
% average_firing = squeeze(mean(spikes_tot_nor));
% std_firing = squeeze(std(spikes_tot_nor)) / sqrt(size(spikes_tot_nor, 1));
average_firing = squeeze(mean(spikes_tot));
std_firing = squeeze(std(spikes_tot)) / sqrt(size(spikes_tot, 1));
subplot(1,2,2)
boundedline(linspace(0, 9, 9000), smoothdata(average_firing, 'gaussian', 500), ...
    smoothdata(std_firing, 'gaussian', 500)) % plot smoothed firing rate
ylabel('firing rate'); xlabel('time (ms)'); hold on
plot([3 3], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
plot([6 6], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')
xlim([0 9])

% %% sort of raster map like github mouseland
% greens = cbrewer('seq', 'Greys', 100); % set fancy colormap
% % spikes_reduced = squeeze(mean(reshape(spikes_tot, size(spikes_tot, 1), 5, []), 2)); % reduce to 50-ms window
% spikes_reduced = squeeze(mean(reshape(spikes_tot, size(spikes_tot, 1), 20, []), 2)); % reduce to 50-ms window
% % sort spike trains by similarity but only use "central" part to stress opto differences
% part_to_use = size(spikes_reduced, 2) / 9 * 2 : size(spikes_reduced, 2) / 9 * 6;
% idx_sorted = sort_spike_trains(zscore(spikes_reduced(:, part_to_use)));
% zscored_units = zscore(spikes_reduced, [], 2); % zscore
% zscored_units(zscored_units > 3) = 3; % take away extreme values for better plotting
% zscored_units(zscored_units < 0) = 0; % take away extreme values for better plotting
% % to have the darkest spike trains during opto on top
% flip_figure = 0;
% if sum(zscored_units(idx_sorted == 1, part_to_use)) > ...
%         sum(zscored_units(idx_sorted == max(idx_sorted), part_to_use))
%     figure; imagesc(zscored_units(idx_sorted, :)); colormap(greens) % plot
% else
%     figure; imagesc(flipud(zscored_units(idx_sorted, :))); colormap(greens) % plot
%     flip_figure = 1; % to have the plot in the same direction also for variance
% end
% hold on
% plot([size(spikes_reduced, 2) / 9 * 3 size(spikes_reduced, 2) / 9 * 3], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
% plot([size(spikes_reduced, 2) / 9 * 6 size(spikes_reduced, 2) / 9 * 6], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
% ylabel('Single units'); xlabel('Time (0.1*s)')

%% sort of raster map like github mouseland
greens = cbrewer('seq', 'Greys', 100); % set fancy colormap
% spikes_reduced = squeeze(mean(reshape(spikes_tot, size(spikes_tot, 1), 5, []), 2)); % reduce to 50-ms window
spikes_reduced = squeeze(mean(reshape(spikes_tot, size(spikes_tot, 1), 20, []), 2)); % reduce to 50-ms window
part_to_use=55:65;
[U, Sv, ~] = svd(zscore(spikes_reduced(:, part_to_use))); % U has "eigenvectors" ("principal components") of S, Sv the "eigenvalues"
[~, idx_sorted] = sort(U(:, 1), 'descend'); % sort by "first component"

figure; imagesc(flipud(zscored_units(idx_sorted, :))); colormap(greens) % plot
flip_figure = 1; % to have the plot in the same direction also for variance


set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')



figure
subplot(1,2,1)
scatter(firing_tot(:, 1), firing_tot(:, 2), 30, 's', 'filled', 'k')
p = signrank(firing_tot(:, 1), firing_tot(:, 2));
hold on
scatter(median(firing_tot(:, 1)), median(firing_tot(:, 2)), 100, 's', 'filled', 'r')
title(strcat('firng rate pre-during p=',num2str(p)))
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca, 'FontName', 'Arial')
xl = xlim; yl = ylim; xlim([min([xl(1) yl(1)]) max([xl(2) yl(2)])]); ylim(xlim);
refline(1,0)
subplot(1,2,2)
scatter(firing_tot(:, 2), firing_tot(:, 3), 30, 's', 'filled', 'k')
p = signrank(firing_tot(:, 2), firing_tot(:, 3));
hold on
scatter(median(firing_tot(:, 2)), median(firing_tot(:, 3)), 100, 's', 'filled', 'r')
title(strcat('firng rate during-post p=',num2str(p)))
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca, 'FontName', 'Arial')
xl = xlim; yl = ylim; xlim([min([xl(1) yl(1)]) max([xl(2) yl(2)])]); ylim(xlim);
refline(1,0)

if strcmp(expType, 'VIR')
    modulated_channels(1, 1) = nnz(pvalue < 0.01 & (stim_layer == 1)' & channels < 25 ...
        & OMI < 0) ./ nnz((stim_layer == 1)' & channels < 25 & ...
        (OMI < 0 | (OMI > 0 & pvalue > 0.01))) * 100;
    modulated_channels(1, 2) = nnz(pvalue < 0.01 & (stim_layer == 1)' & channels > 24 ...
        & OMI < 0) ./ nnz((stim_layer == 1)' & channels > 24 & ...
        (OMI < 0 | (OMI > 0 & pvalue > 0.01))) * 100;
    modulated_channels(2, 2) = nnz(pvalue < 0.01 & (stim_layer == 2)' & channels > 24 ...
        & OMI < 0) ./ nnz((stim_layer == 2)' & channels > 24 & ...
        (OMI < 0 | (OMI > 0 & pvalue > 0.01))) * 100;
    modulated_channels(2, 1) = nnz(pvalue < 0.01 & (stim_layer == 2)' & channels < 25 ...
        & OMI < 0) ./ nnz((stim_layer == 2)' & channels < 25 & ...
        (OMI < 0 | (OMI > 0 & pvalue > 0.01))) * 100;
    modulated_channels
    figure; imagesc(modulated_channels, [10 60]); colormap(greens)
end


display(strcat('there are a total of_', num2str(nnz(pvalue < 0.01 ...
    & OMI > 0)), '_units with pvalue < 0.01 and OMI > 0 out of_', ...
    num2str(length(OMI))))
display(strcat('there are a total of_', num2str(nnz(pvalue < 0.01 ...
    & OMI < 0)), '_units with pvalue < 0.01 and OMI < 0 out of_', ...
    num2str(length(OMI))))
end
