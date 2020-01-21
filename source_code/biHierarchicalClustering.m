%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  biHierarchicalClustering function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs top-down hierarchical clustering in two levels %
% as described in section 'Model-free higher-order drug combination     %
% analysis' and in particular subsection 'top-down hierarchical         %
% clustering' of the 'Results' in the main article text.                %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% r_Q: matrix [N_Txd], where N_T is the total number of treatments and  %
% d is the total number of measured protein releases (after QC). Each   %
% element r_Q(i, j) corresponds to the ratio defined in question Q2 or  %
% Q3.                                                                   %                               
%                                                                       %
% N_W: number of treatments (i.e., exprimental wells) to be clustered.  %
%                                                                       %
% annot_W: cell array {N_Wx1}. A particular cell {i} contains the anno- %
% tation for the corresponding treatment/cell state i.                  %
%                                                                       %
% R: number of repetitions/restarts for K-Means clustering (see 'main.m'%
% where R is set to 30, this can be changed).                           %
%                                                                       %
% tau_SSE_drop: drop in the sum of squared errors (SSE) when transitio- %
% ning from k-1 to k as the number of clusters to be used (see function %
% 'KMeansOptimal' for more details).                                    %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% r_Clusters_annot: cell array {Nx2}, where N is the total number of    %
% different clusters in the 1st hierarchical level. A particular cell   %
% {i, 1} contains a cell array {Mx1}, where M is the total number of    %
% treatments belonging to cluster i. A particular cell array {i, 2}     %
% contains a nested {1xK} cell array, where K is the total number of    %
% subclusters for cluster i. For each subcluster k the corresponding    %
% treatment annotations are included in a {N_T_Sx1} cell array, where   %
% N_T_S is the number of treatments in subcluster k of cluster i.       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==============================================================================================================================================

function r_Clusters_annot = biHierarchicalClustering(r_Q, N_W, annot_W, R, tau_SSE_drop)

%% Perform 2-Level Hierarchical Clustering
[~, ~, ~, idx_opt, k_opt] = KMeansOptimal(r_Q, N_W, R, tau_SSE_drop);
idx_sub_opt = cell(k_opt, 1);
k_sub_opt = zeros(k_opt, 1);
r_Clusters = cell(k_opt, 2);
r_Clusters_annot = cell(k_opt, 2);
c = 0;
for i = 1 : k_opt
   
    c = c + 1;
    
    ind_cluster = (idx_opt == i); % find all indices per cluster
    r_Clusters{i, 1} = r_Q(ind_cluster, :);
    r_Clusters_annot{i, 1} = annot_W(ind_cluster);
    
    if (size(r_Clusters_annot{i, 1}, 1) > 2)  % Perform clustering in the 2nd level if the corresponding 1st cluster contains more than 2 elements
        
        
        [~, ~, ~, idx_sub_opt{i}, k_sub_opt(i)] = KMeansOptimal(r_Q(ind_cluster, :), sum(ind_cluster), R, -0.3);
        for j = 1 : k_sub_opt(i)

            ind_sub_cluster = (idx_sub_opt{i} == j);
            r_Clusters{i, 2}{j} = r_Clusters{i, 1}(ind_sub_cluster, :);
            r_Clusters_annot{i, 2}{j} = r_Clusters_annot{i, 1}(ind_sub_cluster, :);

        end
        
    end
     
end
    
end

%==============================================================================================================================================


