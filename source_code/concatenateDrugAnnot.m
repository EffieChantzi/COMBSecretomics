
%%  concatenateDrugAnnot function - 19/10/29  %%

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%                 Chantzi Effie                 %%
           %%                COMBSecretomics                %%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function merges the names of all drugs of a particular cluster   %
% in one unique string that is used for the corresponding legend in the %
% hierarchical clustering figure after validation.                      %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% A: cell array {Nx1}, where N is the total number of treatments that   %
% belong to a particular cluster. Each cell {i} contains the name of    %
% treatment i in the form 'TX'. Here, 'T' stands for treament and 'X' is%
% the numeric identifier for 'T' (for details see section 'Example raw  %
% data file', Supplement). For example, let us assume that one cluster  %
% consists of the two single treatments T1, T2 and the combination T123.%
% Then the input A would be: A{1} = 'T1', A{2} = 'T2' and A{3} = 'T123'.%
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% A_str: comma-separated string with all drug names per cluster. Using  %
% the aforementioned example, A_str = 'T1,T2,T123'.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===============================================================================

function A_str = concatenateDrugAnnot(A)

N = length(A);

A_str = '';
for i = 1 : N
    
    if (i == N)
        
        A_str = strcat(A_str, A{i});
        
    else

        A_str = strcat(A_str, A{i}, ',');
        
    end
    
end

end

%===============================================================================

