
%%  HierarchicalPartitionResampling function - 19/10/29  %%

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%                 Chantzi Effie                 %%
           %%                COMBSecretomics                %%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds and plots the most common unique partitions in the%
% 1st and 2nd hierarchical levels among all resampling-based validation %
% datasets. For more details, see section 'Resampling statistics' of the%
% 'Results' in the main article.                                        %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% r_Clusters_annot_resampl: {1xN_resampl} cell array containing the     %
% treatment annotations of all (sub)custers in the 1st/2nd hierarchical %
% levels for all N_resampl validation datasets. Each cell {i} contains a%
% nested cell array {N_ix2}, where N_i is the number of clusters in the %
% 1st hierarchical level for validation dataset i. A particular nested  %
% cell array {i}{j, 1} contains the treatment annotations for cluster j %
% as a cell array {N_Tx1}, where N_T is the number of treatments in j. A%
% particular nested cell array {i}{j, 2} contains a cell array {1xK},   %
% where K is the number of subclusters for cluster j. For each sub-     %
% cluster k the corresponding treatment annotations are included in a   %
% {N_T_Sx1} cell array, where N_T_S is the number of treatments in      %
% subcluster k of cluster j. If there is only one dataset, then the     %
% input r_Clusters_annot is a {Nx2} cell array with exactly the same    %
% structure as {N_ix2}.                                                 %
%                                                                       %
% N_resampl: user-defined number of validation datasets to be created   %
% based on resampling.                                                  %
%                                                                       %
% extra_str_title: extra string to be added in the title of the figure  %
% to be generated. This is the number of resamplings to be performed    %
% and/or the identifier for the stimulation used (S1, S2, S3 etc).      %
% However, this can be modified according to the users' preferences.    %
%                                                                       %
% filename_str: filename for saving the figure to be produced.          %
%                                                                       %
% resDir: directory where the generated results should be saved.        %
%                                                                       %
% codeDir: directory with the source code.                              %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% r_Clusters_annot_resampl: {1xN_resampl_unique} cell array containing  %
% the treatment annotations of all unique(sub)custers in the 1st/2nd    %
% hierarchical levels among all rasampling based validation datasets.   %
% Each cell {i} contains a nested cell array {N_ix2}, where N_i is the  %
% number of clusters in the 1st hierarchical level for validation       %
% dataset i. A particular nested cell array {i}{j, 1} contains the      %
% treatment annotations for cluster j as a cell array {N_Tx1}, where N_T%
% denotes the number of treatments in j. A particular nested cell array %
% {i}{j, 2} contains a cell array {1xK} where K is the number of        %
% subclusters for cluster j. For each subcluster k the corresponding    %
% treatment annotations are included in a {N_T_Sx1} cell array, where   %
% N_T_S is the number o treatments in subcluster k of cluster j. Notably%
% all cell arrays of r_Clusters_annot_resampl are equivalent in content %
% meaning that the corresponding partitions are identical. However, the %
% order of these partitions might differ. For example, let us assume    %
% that the 1st cell array contains the following two 1st level clusters:%
% r_Clusters_annot_resampl{1}{1, 1}(:, 1)= {'T23'}                      %
% r_Clusters_annot_resampl{1}{2, 1}(:, 1)= {'T1','T13','T2','T3','T12'} %
% while the 2nd cell array contains:                                    %
% r_Clusters_annot_resampl{2}{1, 1}(:, 1)= {'T1','T13','T2','T3','T12'} % 
% r_Clusters_annot_resampl{2}{2, 1}(:, 1)= {'T23'}.                     %
% The partitions demonstrated by r_Clusters_annot_resampl{1} and        %
% r_Clusters_annot_resampl{2} are equivalent in content.                %
%                                                                       %
% Figure showing the normalized frequency (%) of all unique partitions  %
% in the 1st hierarchical level among all resampling based validation   %
% datasets. Annotations are shown in the figure only for the 3  most    %
% frequent partitions. The figure is saved in the directory resDir (see %
% above).                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=================================================================================================================================================

function r_Clusters_annot_resampl = HierarchicalPartitionResampling(r_Clusters_annot_resampl, N_resampl, extra_str_title, filename_str, resDir, ...
                                                                                                                                        codeDir)

[r_Clusters_annot_resampl_simplified_1, r_Clusters_annot_resampl_simplified_2] = concatenateClusterAnnot(r_Clusters_annot_resampl);
[r_Clusters_annot_resampl_unique, r_Clusters_annot_resampl_occur] = findClusterOccurence(r_Clusters_annot_resampl_simplified_1, '', []);
r_Clusters_annot_resampl_occur = 100*(r_Clusters_annot_resampl_occur/N_resampl);
[sort_annot_resampl_occur, ind_sort] = sort(r_Clusters_annot_resampl_occur, 'descend');

load 'cool_colormap.mat' 'cool_colormap';
cmap = cool_colormap;
figure();
plot(r_Clusters_annot_resampl_occur, 'p', 'MarkerEdgeColor', cmap(60, :), 'MarkerFaceColor', cmap(60, :), 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
grid on;
ylim([0 100]);
if (isempty(extra_str_title))
    title_str = sprintf('Hierarchical K-means clustering (1^{st} level, %d resamplings)', N_resampl);
else
    title_str = sprintf('Hierarchical K-means clustering (1^{st} level, %d resamplings, %s)', N_resampl, extra_str_title);
end
title(title_str, 'FontWeight', 'Bold', 'FontSize', 12, 'FontName', 'Sans Serif');
xlabel('unique partitions', 'FontWeight', 'Bold', 'FontSize', 12, 'FontName', 'Sans Serif');
ylabel('Occurence (%)', 'FontWeight', 'Bold', 'FontSize', 12, 'FontName', 'Sans Serif');
for i = 1 : 3
    text_str = r_Clusters_annot_resampl_unique(ind_sort(i));
    text(ind_sort(i), sort_annot_resampl_occur(i) + 5, text_str{:}, 'FontSize', 8, 'FontWeight', 'Bold', 'FontName', 'Sans Serif');
    if (i == 1)
           plot(ind_sort(i), sort_annot_resampl_occur(i), 'p', 'MarkerEdgeColor', cmap(60, :), 'MarkerFaceColor', cmap(3, :), ...
                                                                                                        'LineWidth', 1, 'MarkerSize', 12);
    end
end
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat(filename_str, '.eps'));
clear tmp;
cd(codeDir);
  
%% Most common partition in the 1st Hierarchical Level
winner_1 = r_Clusters_annot_resampl_unique(ind_sort(1));
winner_1 = winner_1{:};
ind_pivot = zeros(1, N_resampl);
for k = 1 : N_resampl
    intersection = intersect(r_Clusters_annot_resampl_simplified_1{k}, winner_1);
    if (length(intersection) == length(winner_1))

        ind_pivot(k) = 1;

    end
end
ind_pivot = logical(ind_pivot);
r_Clusters_annot_resampl = r_Clusters_annot_resampl(ind_pivot);

%% Unique partitions in the 2nd Hierarchical Level of the most unique partitions in the 1st Hierarchical Level
r_Clusters_annot_resampl_simplified_2 = r_Clusters_annot_resampl_simplified_2(ind_pivot);
ind_sub = cell2mat(cellfun(@(x) size(x, 2), r_Clusters_annot_resampl_simplified_2, 'UniformOutput', false));
ind_sub_unique = unique(ind_sub, 'stable');
no_ind_sub_unique = length(ind_sub_unique);
sum_count = zeros(1, no_ind_sub_unique);
for u = 1 : no_ind_sub_unique
    sum_count(u) = sum(ind_sub == ind_sub_unique(u));
end
[sum_count_max, ind_sum_count_max] = max(sum_count);
r_Clusters_annot_resampl_simplified_2 = r_Clusters_annot_resampl_simplified_2(ind_sub == ind_sub_unique(ind_sum_count_max));
r_Clusters_annot_resampl_simplified_2_mod = cell(size(r_Clusters_annot_resampl_simplified_2));
for u = 1 : sum_count_max

    tmp = r_Clusters_annot_resampl_simplified_2{u};
    reshape_dim = cell2mat(cellfun(@(x) size(x, 1), tmp, 'UniformOutput', false));
    tmp = [tmp{:}];
    tmp = reshape(tmp, [sum(reshape_dim) 1]);
    r_Clusters_annot_resampl_simplified_2_mod{u} = tmp;

end

%% Most common partition in the 2nd Hierarchical Level
[r_Clusters_annot_resampl_unique, r_Clusters_annot_resampl_occur] = findClusterOccurence(r_Clusters_annot_resampl_simplified_2_mod, '', []);
[~, ind_max_occur] = max(r_Clusters_annot_resampl_occur);
winner_2 = r_Clusters_annot_resampl_unique(ind_max_occur(1));
winner_2 = winner_2{:};
ind_pivot = zeros(1, sum_count_max);
for u = 1 : sum_count_max
    intersection = intersect(r_Clusters_annot_resampl_simplified_2_mod{u}, winner_2);
    if (length(intersection) == length(winner_2))

        ind_pivot(u) = 1;

    end
end
ind_pivot = logical(ind_pivot);

%% Most common partitions in the 1st and 2nd Hierarchical Levels
r_Clusters_annot_resampl = r_Clusters_annot_resampl(ind_pivot);

end

%=================================================================================================================================================

