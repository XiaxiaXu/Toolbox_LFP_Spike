function zscored_units=spike_sort_order(SpikesMatrix,part_to_use)
% Xiaxia
SpikesMatrix_toUse=SpikesMatrix(:, part_to_use); 
[U, Sv, ~] = svd((SpikesMatrix_toUse)); % U has "eigenvectors" ("principal components") of S, Sv the "eigenvalues"
[~, idx_sorted] = sort(U(:, 1), 'descend'); % sort by "first component"
zscored_units = zscore(SpikesMatrix, [], 2); % zscore
zscored_units=zscored_units(idx_sorted, :);

