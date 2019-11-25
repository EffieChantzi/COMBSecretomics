%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2019  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  exhaustiveSubsetSearch function - 19/10/29  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the exhaustive subset search for the top-down  %
% hierarchical clustering as described in the subsection 'Top-down      %
% hierarchical clustering' of the 'Results' in the main article.        %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% l: cell array {Nx1}, where N is the number of treatments belonging to %
% a particular (sub)cluster. A particular cell {i} contains the numeric %
% part of the annotation of treatment i (for details see section        %
% 'Example raw data file', Supplement). For example, let us assume that %
% a cluster consists of the single treatment T1 and the two pairs T23   %
% T13. Then l would be: l{1} = '1', l{2} = '23' and l{3} = '13'.        %
%                                                                       %
% l_unique: cell array with the unique treatments per (sub)cluster after%     
% each iteration.                                                       %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% l_unique: cell array {Mx1}, where M is the number of unique treatments%
% per (sub)cluster after the exhaustive subset search is terminated. A  %
% particular cell {i} contains the numeric identifier for unique        %
% treatment i.                                                          %
%                                                                       %
% Using the previous example, the corresponding function call would be: %
% >>  l_unique = exhaustiveSubsetSearch(l, [])                          %
%                                                                       %
%  with output l_unique being:                                          %
% l_unique{1} = '1', l_unique{2} = '23' .                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%============================================================================================================================================

function l_unique = exhaustiveSubsetSearch(l, l_unique)

if (isempty(l))
    
    return;
    
else
    
    char_length_l = cell2mat(cellfun(@(x) length(x), l, 'UniformOutput', false));
    ind_l_pivot = (char_length_l == min(char_length_l));
    l_pivot = l(ind_l_pivot);
    l_for_pivot = l(~ind_l_pivot);
    N_l_pivot = length(l_pivot);

    l_unique = [l_unique l_pivot];
    if (isempty(l_for_pivot))

        return;

    else

        for i = 1 : N_l_pivot

            ind_contains = cell2mat(cellfun(@(x, y) sum(ismember(x, y)) > 0, l_for_pivot, repmat(cellstr(l_pivot{i}), size(l_for_pivot)),  ...
                                                                                                                    'UniformOutput', false));
            l_for_pivot(ind_contains) = [];

        end
        l_unique = exhaustiveSubsetSearch(l_for_pivot, l_unique);

    end


end

%============================================================================================================================================

