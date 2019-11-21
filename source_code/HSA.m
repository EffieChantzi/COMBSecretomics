
%% HSA function - 19/10/29  %%

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%                 Chantzi Effie                 %%
           %%                COMBSecretomics                %%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the generalized highest single agent analysis  %
% (HSA) in terms of multi-dimensional protein release patterns for any  %
% combination order. See section 'Generalization of the highest single  %
% agent principle' of the Results in the main article for more details. %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% r_Q: matrix [N_Txd], where N_T is the total number of treatments and  %
% d is the total number of measured protein releases (after QC). Each   %
% element r_Q(i, j) corresponds to the ratio defined in question Q2 or  %
% Q3.                                                                   % 
%                                                                       %
% N_T: number of treatments/experimental wells (i.e., rows of r_Q).     %
%                                                                       %
% annot_W: cell array with as many cells as N_T (see above). Each cell  %
% cell {i} contains the annotation for the corresponding treatment of   %
% in well i. For more details, see section 'Example raw data file' of   %
% the Supplement.                                                       %
%                                                                       %
% extra_str_title: extra string to be added in the title of the HSA     %
% figure. This is the identifier for the stimulation used (S1, S2, S3   %
% etc). In this way, one can disentangle the figures obtained for the   %
% different stimulations. However, this input can be changed according  %
% to the users' preferences.                                            %
%                                                                       %
% question_ID: string that is attached to the filename of the figure    %
% when saved. This could be 'Q2', or 'Q3' based on the first input r_Q  %
% (see above). However, this input can be changed according to the      %
% users' preferences.                                                   %
%                                                                       %
% resDir: directory where the generated results should be saved.        %
%                                                                       %
% codeDir: directory with the source code.                              %
%                                                                       %
% flag_fig: numeric identifier, which must be set to 1 for generating   %
% the HSA figure and 0 otherwise.                                       %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% HSA_r_Q: HSA indices for all combination treatments included in the   %
% experiment as a [N_higher_orderx1] column vector. N_higher_order is   %
% the number of different combination treatments used (see below).      %
%                                                                       %
% annot_W_HSA: annotations for combination treatments as a cell array   %
% {N_higher_orderx1}. A particular cell {i} contains the corresponding  %
% annotation for the combination treatment in well i. For more details, %
% see section 'Example raw data file' in the Supplementary Information. %
%                                                                       %
% N_higher_order: rows of output HSA_r_Q (see above), meaning the number%
% of combination treatments used.                                       %
%                                                                       %
% !! If fig_flag set to 1, a heatmap with the HSA indices for all       %
% N_higher_order combination treatments is generated and saved under    %
% resDir (see above).                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==============================================================================================================================================

function [HSA_r_Q, annot_W_HSA, N_higher_order] = HSA(r_Q, N_T, annot_W, extra_str_title, question_ID, resDir, codeDir, fig_flag)

if (fig_flag ~= 0) && (fig_flag ~= 1)
   
    error('Invalid flag for figure.');
    
end

r_Q = abs(r_Q);
N_P = size(r_Q, 2);
r_Q = sum(r_Q, 2)/N_P;

annot_W_HSA = cellfun(@(x) char(x), annot_W, 'UniformOutput', false);
annot_W_HSA = cellfun(@(x) x(2 : end), annot_W_HSA, 'UniformOutput', false);
    
first_order_ind = cell2mat(cellfun(@(x) length(x) == 1, annot_W_HSA, 'UniformOutput', false));
N_first_order = sum(first_order_ind);
N_higher_order = N_T - N_first_order;

higher_order_annot = annot_W_HSA(~first_order_ind);
higher_order_annot = cell2mat(cellfun(@(x) str2double(x), higher_order_annot, 'UniformOutput', false));
[~, I] = sort(higher_order_annot, 'ascend');
higher_order_annot = higher_order_annot(I);

HSA_r_Q = zeros(N_higher_order, 1);
annot_W_HSA = cell(1, N_higher_order);
for i = 1 : N_higher_order

    higher_order_annot_tmp = num2str(higher_order_annot(i));
    s = extractSubsets(higher_order_annot_tmp);
    s = cellfun(@(x) strcat('T', x), s, 'UniformOutput', false);
    N_subsets = size(s, 1);
    r_Q_for_min = zeros(N_subsets, 1);
    for j = 1 : N_subsets

        ind = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W, repmat(cellstr(s{j}), N_T, 1), 'UniformOutput', false));
        r_Q_for_min(j) = r_Q(ind);

    end
    r_Q_min = min(r_Q_for_min); 
    
    higher_order_annot_tmp = strcat('T', num2str(higher_order_annot(i)));
    ind_higher_order = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W, repmat(cellstr(higher_order_annot_tmp), N_T, 1), ...
                                                                                                        'UniformOutput', false));
    HSA_r_Q(i) = r_Q_min - r_Q(ind_higher_order);
    annot_W_HSA{i} = higher_order_annot_tmp;

end

if (fig_flag)
    
    figure();
    h = heatmap({'all proteins'}, annot_W_HSA, round(HSA_r_Q, 2, 'significant'));
    if (isempty(extra_str_title))
        h.Title = 'HSA Index';
    else
        h.Title = sprintf('HSA Index (%s)', extra_str_title);   
    end
    h.Colormap = cool;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    cd(resDir);
    print -depsc -painters -r400 tmp
    if (isempty(extra_str_title))
        filename = strcat('HSA_', question_ID, '.eps');
    else
        filename = strcat('HSA_', extra_str_title, '_', question_ID, '.eps');
    end
    movefile('tmp.eps', filename);
    clear tmp;
    cd(codeDir);

end

end

%==============================================================================================================================================


