%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  findClusterOccurence function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the unique partitions in each hierarchical level  %
% as obtained from the top-down hierarchical K-Means clustering for all %
% different resampling-based validation datasets. For more details, see %
% section 'Resampling statistics' of 'Results' in the main article text.%
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% l: {1xN_resampl} cell array with the partitions of the 1st or 2nd     %
% hierarchical levels for all N_resampl validation datasets. Each cell  %
% {i} contains a {N_ix1} cell array with each element {i}{j} being a    %
% concatenated comma separated string of all treatments in (sub)cluster %
% j for validation dataset i.                                           %
%                                                                       %
% l_unique: {1xN_resampl_unique_cur} cell array with the unique         %
% partitions of the 1st or 2nd hierarchical levels after each iteration.%
% (i.e., currently). Each cell {i} contains a {N_ix1} cell array with   %
% each element {i}{j} being a concatenated comma separated string of all%
% treatments in the current unique (sub)cluster j.                      %                                                  
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% l_unique: {1xN_resampl_unique} cell array with the unique partitions  %
% of the 1st or 2nd hierarchical levels after all iterations at the end %
% of the search. Each cell {i} contains a {N_ix1} cell array with each  %
% element {i}{j} being a concatenated comma separated string of all     %
% treatments in the unique (sub)cluster j.                              %
%                                                                       %
% l_unique_occur: row vector [1xN_resampl_unique], where each element   %
% (i) contains the number of occurrences for the unique partition i of  %
% the output l_unique across input l (see above).                       %
%                                                                       %
% 'findClusterOccurence' is a recursive type of function and should be  %
% initiated by the following call:                                      %
% >> [l_unique, l_unique_occur] = findClusterOccurence(l, '', []);      %
%                                                                       %
% See function 'HierarchicalPartitionResampling' for more details on    %
% how 'findClusterOccurence' is employed by COMBSecretomics.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%============================================================================================================

function [l_unique, l_unique_occur] = findClusterOccurence(l, l_unique, l_unique_occur)

if (isempty(l))
   
    return;
    
else
    
    l_pivot = l{1};
    l_for_pivot = l(2 : end);
    N_l_for_pivot = length(l_for_pivot);
    l_pivot_cell{1} = l_pivot;
    
    if (N_l_for_pivot == 1)
        
        l_unique = [l_unique l_pivot_cell];
        occur = 1;
        l_unique_occur = [l_unique_occur occur];
        l_for_pivot = [];
    
    else
       
        ind_pivot = zeros(1, N_l_for_pivot);
        for i = 1 : N_l_for_pivot
            intersection = intersect(l_for_pivot{i}, l_pivot);
            if (length(intersection) == length(l_pivot))

                ind_pivot(i) = 1;

            end

        end
        ind_pivot = logical(ind_pivot);
        l_unique = [l_unique l_pivot_cell];
        occur = sum(ind_pivot) + 1;
        l_unique_occur = [l_unique_occur occur];
        l_for_pivot = l_for_pivot(~ind_pivot);
        
    end
    
    [l_unique, l_unique_occur] = findClusterOccurence(l_for_pivot, l_unique, l_unique_occur);
    
end

end

%============================================================================================================

