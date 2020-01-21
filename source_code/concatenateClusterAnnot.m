%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  concatenateClusterAnnot function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function concatenates the treatment names per (sub)cluster as one%
% unified comma separated string.                                       %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% r_Clusters_annot: {1xN_resampl} cell array containing the treatment   %
% annotations of all (sub)custers in the 1st and 2nd hierarchical levels%
% for all N_resampl validation datasets. Each cell {i} contains a nested%
% cell array {N_ix2}, where N_i is the number of clusters in the 1st    %
% hierarchical level for validation dataset i. A particular nested cell %
% array {i}{j, 1} contains the treatment annotations for cluster j as a %
% cell array {N_Tx1}, where N_T is the number of treatments in j. A     %
% particular nested cell array {i}{j, 2} contains a cell array {1xK},   %
% where K is the number of subclusters for cluster j. For each sub-     %
% cluster k the corresponding treatment annotations are included in a   %
% {N_T_Sx1} cell array, where N_T_S is the number of treatments in      %
% subcluster k of cluster j. If there is only one dataset, then the     %
% input r_Clusters_annot is a {Nx2} cell array with exactly the same    %
% structure as {N_ix2}.                                                 %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% r_Clusters_annot_simplified_1: {1xN_resampl} cell array with the 1st  %
% level partitions for all N_resampl validation datasets. Each cell {i} %
% contains a nested cell array {N_ix1} with each element {i}{j} being a %
% concatenated comma separated string of all treatments in cluster j in %
% the 1st hierarchical level for validation dataset i. If there is only %
% one dataset, then the output r_Clusters_annot_simplified_1 is a {Nx1} %
% cell array with the same structure as {N_ix1}.                        % 
%                                                                       %
% r_Clusters_annot_simplified_2: {1xN_resampl} cell array with the 2nd  %
% level partitions all N_resampl validation datasets. Each cell {i}     %
% contains a nested cell array {1xK_i}, where K_i denotes the number of %
% clusters in the 1st hierarchical level for validation dataset i. Each %
% element {i}{j} is a {N_T_Sx1} cell array with N_T_S being the number  %
% of subclusters for cluster j and each element {i}{j}{l} being a       %
% concatenated comma separated string of all treatments in cluster j    %
% that belong to the corresponding subcluster l. If there is only one   %
% dataset, then the output r_Clusters_annot_simplified_2 is a {1xK} cell%
% array with the same structure as {1xK_i}.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%====================================================================================================================

function [r_Clusters_annot_simplified_1, r_Clusters_annot_simplified_2] = concatenateClusterAnnot(r_Clusters_annot)

N_rows = size(r_Clusters_annot, 1);
if (N_rows == 1)
    
    N = size(r_Clusters_annot, 2);
    r_Clusters_annot_simplified_1 = cell(size(r_Clusters_annot));
    r_Clusters_annot_simplified_2 = cell(size(r_Clusters_annot));
    
else
    
    N = 1;
    r_Clusters_annot{1} = r_Clusters_annot(:, 1);

end

r_Clusters_annot_1 = cellfun(@(x) x(:, 1), r_Clusters_annot, 'UniformOutput', false); 
r_Clusters_annot_2 = cellfun(@(x) x(:, 2), r_Clusters_annot, 'UniformOutput', false); 
for i = 1 : N
   
    tmp_1 = r_Clusters_annot_1{i};
    tmp_1 = cellfun(@(x) [x{:}], tmp_1, 'UniformOutput', false);
    ind = cell2mat(cellfun(@(x) iscell(x), tmp_1, 'UniformOutput', false));
    if (sum(ind) == 0)
        tmp_1 = r_Clusters_annot_1{i};
    end
    tmp_1 = cellfun(@(x) strjoin(x, ','), tmp_1, 'UniformOutput', false);
   
    if (N_rows == 1)
        
        r_Clusters_annot_simplified_1{i} = tmp_1;
        tmp_2 = r_Clusters_annot_2{i};
        tmp_2 = cellfun(@(x) x', tmp_2, 'UniformOutput', false);
        no_sub = size(tmp_2, 1);
        tmp_2_tmp = cell(1, no_sub);
        for j = 1 : no_sub
            if (isempty(tmp_2{j}))
               tmp_2_tmp{j} = [];
            else
                y = cellfun(@(x) [x{:}], tmp_2{j}, 'UniformOutput', false);
                ind_y = cell2mat(cellfun(@(x) iscell(x), y, 'UniformOutput', false));
                if (sum(ind_y) == 0)
                    y = tmp_2{j};
                end
                tmp_2_tmp{j} = cellfun(@(x) strjoin(x, ','), y, 'UniformOutput', false);
            end
        end
        r_Clusters_annot_simplified_2{i} = tmp_2_tmp;
        
    else
        
        r_Clusters_annot_simplified_1 = tmp_1;
        r_Clusters_annot_simplified_2 = [];
        
    end
    
end

end

%====================================================================================================================

