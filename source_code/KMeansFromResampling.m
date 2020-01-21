%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  HierarchicalPartitionResampling function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function retrieves the protein release patterns for the most     %
% common partitions in the 1st/2nd hierarchical levels (for more details%
% see function 'HierarchicalPartitionResampling'), which are then used  %
% by the function 'plotHierarchicalKMeans' to calculate and plot the    %
% corresponding centroids.                                              %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% r_Clusters_annot_winner: cell array {Nx2} containing the annotations  %
% for the most frequent partition in the 1st/2nd hierarchical levels. N %
% where N denotes the number of clusters in the 1st hierarchical level. %
% Each cell {i, 1} contains a cell array {Mx1} with M being the total   %
% number of treatments belonging to cluster i. A particular cell array  %
% {i, 2} contains a nested {1xK} cell array, where K is the total number%
% of subclusters for cluster i. For each subcluster k the corresponding %
% treatment annotations are included in a {N_T_Sx1} cell array, where   %
% N_T_S is the number of treatments in subcluster k of cluster i.       %
%                                                                       %
% r_Q: matrix [N_Txd], where N_T is the total number of treatments and  %
% d is the total number of measured protein releases (after QC). Each   %
% element r_Q(i, j) corresponds to the ratio defined in question Q2 or  %
% Q3.                                                                   %
%                                                                       %
% N_T: number of treatments/experimental wells (i.e., rows of r_Q).     %
%                                                                       %
% N_P: number of measured proteins per treatment/experimental well.     %
%                                                                       %
% annot_W_unique_treatments: cell array with as many cells as the number%
% of unique treatments. A particular cell {i} contains the annotation   %
% for treatment i. For more details, see section 'Example raw data file'%
% of the Supplement.                                                    %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% r_Clusters_winner: {Nx2} cell array containing the values r_Q for the %
% most frequent partition in the 1st/2nd hierarchical levels. Each cell %
% array {i, 1} contains a matrix [N_T_CxN_P] (part of input r_Q), where %
% N_T_C is the number of treatments in cluster i. A particular cell     %
% array {i, 2} contains a nested {1xK} cell array, where K is the total %
% number of subclusters for cluster i. A particular nested cell array   %
% {i, 2}{1, k} contains a matrix [N_T_SxN_P], where N_T_S denotes the   %
% number of treatments in subcluster k of cluster i.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=========================================================================================================================================

function r_Clusters_winner = KMeansFromResampling(r_Clusters_annot_winner, r_Q, N_T, N_P, annot_W_unique_treatments)

r_Clusters_winner = cell(size(r_Clusters_annot_winner));
K = size(r_Clusters_annot_winner, 1);
for i = 1 : K

    for j = 1 : 2

        tmp = r_Clusters_annot_winner{i, j};
        if (j == 1)

            no_data_points = size(tmp, 1);
            val = zeros(no_data_points, N_P);
            for d = 1 : no_data_points
                tmp_d = tmp{d};
                if ischar(tmp_d)
                    tmp_d = cellstr(tmp_d);
                end
                ind = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W_unique_treatments, repmat(tmp_d, N_T, 1), 'UniformOutput', false));
                val(d, :) = r_Q(ind, :);
            end
            r_Clusters_winner{i, j} = val;

        else

            no_sub_clusters = size(tmp, 2);
            for s = 1 : no_sub_clusters

                tmp_sub_cluster = tmp{s};
                no_data_points = size(tmp_sub_cluster, 1);
                val = zeros(no_data_points, N_P);
                for d = 1 : no_data_points
                    tmp_sub_cluster_d = tmp_sub_cluster{d};
                    if ischar(tmp_sub_cluster_d)
                        tmp_sub_cluster_d = cellstr(tmp_sub_cluster_d);
                    end
                    ind = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W_unique_treatments, repmat(tmp_sub_cluster_d, N_T, 1), ...
                                                                                                               'UniformOutput', false));
                    val(d, :) = r_Q(ind, :);
                end
                r_Clusters_winner{i, j}{s} = val;

            end

        end

    end

end

end

%=========================================================================================================================================

