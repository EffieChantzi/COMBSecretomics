%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  extractSubsets function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extracts all possible lower-order subsets of a higher-  %
% order combination treatment. This function is needed for employing    %
% the generalized highest single agent principle according to Eq. (2)   %
% in the main article. The numeric annotation for the treatments must   %
% follow the instructions given in section 'Example raw data file' of   %
% the Supplementary Information.                                        %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% T_name: numeric identifier for treatment (for details see section     %
% 'Example raw data file', Supplementary Information). For example if   %
% one wants to retrieve all possible subsets of the 3-order combination %
% treatment T123, then T_name = '123'.                                  %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% s: cell array {Nx1}, where N is the total number of subsets for the   %
% input treatment T_name. Each cell {i} contains the numeric identifier %
% for subset i. Taking the previous example, the elements of the output %
% s would be:                                                           %
% s{1} = '1', s{2} = '2', s{3} = '3', s{4} = '12', s{5} = '13' and      %
% s{6} = '23'.                                                          %
%                                                                       %
% In order to run the aforementioned example run on MATLAB's command    %
% window:                                                               %
% >> s = extractSubsets('123')                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===============================================================================

function s = extractSubsets(T_name)

no_drugs = length(T_name);

drugs = cell(1, no_drugs);
for i = 1 : no_drugs
   
    drugs{i} = T_name(i);
    
end


count = 0;
for i = 1 : no_drugs - 1
    
    s_tmp = nchoosek(str2double(drugs{1}) : str2double(drugs{end}), i);
    r = size(s_tmp, 1);
    for k = 1 : r
       
        count = count + 1;
        
    end
    
end
N_subsets = count;


s = cell(N_subsets, 1);
count = 0;
for i = 1 : no_drugs - 1
    
    s_tmp = nchoosek(str2double(drugs{1}) : str2double(drugs{end}), i);
    [r, c] = size(s_tmp);
    for k = 1 : r
       
        count = count + 1;
        s_tmp_2 = '';
        for l = 1 : c
        
            s_tmp_2 = strcat(s_tmp_2, num2str(s_tmp(k, l)));
            
        end
        s{count, 1} = s_tmp_2;
        
    end
    
end

end

%===============================================================================

