function spikes_tot = plotSUAfiringrateKlustapulses(Experiment, expType, StimRegion, RespArea,Path)
%Path = get_path;
Gwindow = gausswin(1001, 10); % gaussian window of 1000ms with stdev of 100ms
Gwindow = Gwindow / sum(Gwindow); % normalize the gaussian kernel
save_data = 1;
spikes_tot = [];
firing_tot = [];
OMI = [];
pvalue = [];
pre_stim =30: 40; % in ms, ramp format
stim = 60: 80; % in ms, ramp format
post_stim = 81 : 120; % in ms, ramp format
stim_layer = [];

for n_animal =1: length(Experiment)
    n_animal
    experiment = Experiment(n_animal);
    if nnz( ~ isnan(experiment.StimRegion)) && strcmp(expType,experiment.Exp_type)
        if ~isempty(experiment.animal_ID)
            if strcmp(RespArea, 'PFC2')
                foldersave = strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAPFC23',filesep,experiment.name);
            elseif strcmp(RespArea, 'PFC5')
                foldersave = strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAPFC56',filesep,experiment.name);
            elseif strcmp(RespArea, 'HP')
                foldersave = strcat(Path.output, filesep , 'PulsesStimulationStim',StimRegion,'_SUAHP',filesep,experiment.name);
            end
            
            if exist(strcat(foldersave,'\SUAdata_pulses.mat'))
                load(strcat(foldersave,'\SUAdata_pulses.mat'))
                spikes_matrix = SUAdata_pulses.pulse_spike_matrix;
                if size(spikes_matrix, 2) > 1
                    spikes_units = squeeze(mean(spikes_matrix));
                else
                    spikes_units = squeeze(mean(spikes_matrix))';
                end
                idx=find (sum(spikes_units,2) ==0)
                if sum(idx)~=0
                    spikes_units(idx,:)=[];
                    SUAdata_pulses.OMI(idx)=[];
                    SUAdata_pulses.pvalue(idx)=[];
                end
                firing_units(:, 1) = log10(mean(spikes_units(:, pre_stim), 2));
                firing_units(:, 2) = log10(mean(spikes_units(:, stim), 2));
                firing_units(:, 3) = log10(mean(spikes_units(:, post_stim), 2));
                spikes_tot = cat(1, spikes_tot, spikes_units);
                firing_tot = cat(1, firing_tot, firing_units);
                OMI = horzcat(OMI, SUAdata_pulses.OMI); % concatenate
                pvalue = horzcat(pvalue, SUAdata_pulses.pvalue); % concatenate
                
                clear firing_units
            end
        end
    end
end


figure
RdBu = cbrewer('div', 'RdBu', 100);
num_spikes = sum(spikes_tot, 2); % for size only
num_spikes(num_spikes == 0) = 0.1; % scatter can't handle 0 as an input for size
max_spikes = max(num_spikes);
scaling_factor = 200 / max_spikes;

subplot(2,2,1)
scatter(OMI, - log10(pvalue), num_spikes * scaling_factor, pvalue < 0.001, 'filled');
colormap(flip(RdBu)); caxis([-0.1 1.1])
hold on
plot([0 0], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for no modulation
xlim([-1.1 1.1]); ylim([-0.1 ceil(max(- log10(pvalue)))])
ylabel('- Log10 (pvalue)'); xlabel('Modulation Index stim/pre'); alpha(0.5)
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')

modulation = [nnz(OMI < 0 & pvalue < 0.05) nnz(pvalue > 0.05) nnz(OMI > 0 & pvalue < 0.05)];%% mattia use p<0.05
modulation = modulation ./ length(OMI);

subplot(2,2,2)
bar([modulation; nan(1, length(modulation))], 'stacked')
pie(modulation)
set(gca,'xtick',1); xlim([0 2]);
ylabel('Proportion');% xlabel(''); xticks(''); xticklabels('')
title('Proportion of (un)modulated units stim/pre')
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')

spikes_tot(:,50:55)=[];
spikes_tot=[spikes_tot(:,1:6),spikes_tot];

subplot(2,2,3)
hold on
average_firing = squeeze(mean(spikes_tot));
std_firing = squeeze(std(spikes_tot)) / sqrt(size(spikes_tot, 1));
boundedline(linspace(0, 120, 120), smooth(average_firing, 5),smooth(std_firing, 5)) % plot smoothed firing rate
ylabel('firing rate'); xlabel('time (ms)'); hold on
plot([50 50], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
plot([53 53], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')
xlim([0 120])


subplot(2,2,4)
%spikes_tot_increase=spikes_tot(OMI > 0 & pvalue < 0.05 , :);
spikes_tot_increase=spikes_tot(OMI > 0 , :);
spikes_tot_nor=spikes_tot_increase;
average_firing = squeeze(mean(spikes_tot_nor));

std_firing = squeeze(std(spikes_tot_nor)) / sqrt(size(spikes_tot_nor, 1));
boundedline(linspace(0, 120, 120), smooth(average_firing, 5),smooth(std_firing,  5)) % plot smoothed firing rate
ylabel('firing rate'); xlabel('time (ms)'); hold on
plot([50 50], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
plot([55 55], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')
xlim([0 120])


% greens = cbrewer('seq', 'Greens', 100); % set fancy colormap
% spikes_reduced= squeeze(mean(reshape(spikes_tot, size( spikes_tot, 1), 1, []), 2)); % reduce to 5ms window
% % sort spike trains by similarity but only use "central" part to stress opto differences
% part_to_use = 1:size(spikes_reduced, 2);
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
% plot([size(spikes_reduced, 2) / 2-10  size(spikes_reduced, 2) / 2-10], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
% plot([size(spikes_reduced, 2) / 2-5  size(spikes_reduced, 2) / 2-5], get(gca, 'ylim'), 'k', 'linewidth', 1) % reference line for opto
% ylabel('Single units'); xlabel('Time (0.1*s)')
% set(gca, 'TickDir', 'out'); set(gca, 'FontSize', 14); set(gca, 'FontName', 'Arial')

figure
subplot(1,2,1)
scatter(firing_tot(:, 1), firing_tot(:, 2), 20, 's', 'filled', 'k')
p = signrank(firing_tot(:, 1), firing_tot(:, 2));
hold on
scatter(median(firing_tot(:, 1)), median(firing_tot(:, 2)), 50, 's', 'filled', 'r')
title(strcat('firng rate pre-during p=',num2str(p)))
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca, 'FontName', 'Arial')
xl = xlim; yl = ylim; xlim([min([xl(1) yl(1)]) max([xl(2) yl(2)])]); ylim(xlim);
refline(1,0)
subplot(1,2,2)
scatter(firing_tot(:, 2), firing_tot(:, 3), 20, 's', 'filled', 'k')
p = signrank(firing_tot(:, 2), firing_tot(:, 3));
hold on
scatter(median(firing_tot(:, 2)), median(firing_tot(:, 3)), 50, 's', 'filled', 'r')
title(strcat('firng rate during-post p=',num2str(p)))
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca, 'FontName', 'Arial')
xl = xlim; yl = ylim; xlim([min([xl(1) yl(1)]) max([xl(2) yl(2)])]); ylim(xlim);
refline(1,0)


display(strcat('there are a total of_', num2str(nnz(pvalue < 0.01 ...
    & OMI > 0)), '_units with pvalue < 0.05 and OMI > 0 out of_', ...
    num2str(length(OMI))))
display(strcat('there are a total of_', num2str(nnz(pvalue < 0.05 ...
    & OMI < 0)), '_units with pvalue < 0.05 and OMI < 0 out of_', ...
    num2str(length(OMI))))