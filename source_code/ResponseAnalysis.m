
%%  ResponseAnalysis function - 19/10/29  %%

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%                 Chantzi Effie                 %%
           %%                COMBSecretomics                %%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function normalizes the protein release differences per plate    %
% for stimulated cells. More precisely, this function gives answers to  %
% questions Q1-Q3 described in the 'Results' section of the article     %
% under the subsection 'Normalization of protein release differences'.  %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% F: pre-processed protein release data, meaning after quality control. %
% See function 'preProcessing' for more details. Rows correspond to the % 
% experimental wells (i.e., cell states) and columns to the measured    %
% proteins.                                                             %
%                                                                       %
% annot_W: cell array with as many cells as the number of experimental  %
% wells remained after quality control. A particular cell {i} contains  %
% the annotation for the cell state of the corresponding well i.        %
%                                                                       %
% F_merged: matrix with median values of the raw intra-plate replicate  %
% measurements after quality control. The number of rows correspond to  %
% the number of unique cell states and number of columns to the number  %
% of measured proteins remained after quality control.                  %
%                                                                       %
% annot_W_unique: cell array with as many cells as the number of unique %
% wells after quality control and intra-plate averaging. A particular   %
% cell {i} contains the annotation for the corresponding unique cell    %
% state.                                                                %
%                                                                       %
% N_P: number of measured proteins remained after quality control.      %
%                                                                       %
% annot_P: cell array with as many cells as the number of measured      %
% proteins N_P after quality control. A particular cell {i} contains    %
% the name of the corresponding measured protein i.                     %
%                                                                       %
% N_resampl: user-defined number of resampling-based validation datasets%
% to be created.                                                        %
%                                                                       %
% BarCode: serial barcode of the experimental plate. This is retrieved  %
% automatically from the name of the raw data file, which must be in    %
% the form <barcode.csv>.                                               %
%                                                                       %
% resDir: directory where the generated results should be saved.        %
%                                                                       %
% codeDir: directory with the source code.                              %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% r_Q1: ratio in eq. (1) for answering question Q1 (see main article).  %
% It is stored as a [1xd] row vector where d is the number of measured  %
% protein releases. It defines the total therapeutic need.              %
%                                                                       %
% r_Q2: ratio in eq. (1) for answering question Q2 (see main article).  %
% It is stored as a [Nxd] row vector where N is the number of unique    %
% treatments (i.e., after intra-plate median averaging) and d is the    %
% number of measured protein releases. It defines the modulation        %
% capacities of the tested treatments.                                  %
%                                                                       %
% r_Q3: ratio in eq. (1) for answering question Q3 (see main article).  %
% It is stored as a [Nxd] row vector where N is the number of unique    %
% treatments (i.e., after intra-plate median averaging) and d is the    %
% number of measured protein releases. It defines the restoration       %
% capacities of the tested treatments.                                  %
%                                                                       %
% r_Q1_Q3_S_arranged: matrix containing the ratios r_Q1 (see above) for %
% for all different stimulations. They are stored as a [N_Sxd] matrix,  %
% where N_S is the total number of different stimulations used and d is %
% the number of different proteins measured. As explained in the main   %
% article, these ratios from question Q1 define the therapeutic needs,  %
% which are jointly visualized with the corresponding restoration       %
% capacities defined by r_Q3.                                           %
%                                                                       %
% r_Q3_resampl: cell array {N_SxN_resampl}, where N_S is the total      %
% number of different stimulations used, while N_resampl is the user-   %
% defined number of validation datasets to be created. A particular cell%
% {i, j} contains the matrix r_Q3 (see above) for stimulation i and     %
% validation dataset j.                                                 %
%                                                                       %
% annot_T_S_Q2: cell array {N_Sx1}, where N_S is the total number of    %
% different stimulations used for question Q2. Each cell {i} contains a %
% cell array {N_Tx1}, with the names of all treatments for stimulation  %
% i. N_T is the number of total treatments used for stimulation i.      %
%                                                                       %
% annot_T_S_Q3: cell array {N_Sx1}, where N_S is the total number of    %
% different stimulations used for question Q3. Each cell {i} contains a %
% cell array {N_Tx1}, with the names of all treatments for stimulation  %
% i. N_T is the number of total treatments used for stimulation i.      %
%                                                                       %
% D_str: string identifier for the disease model of interest. This is   %
% automatically retrieved from the raw data file (see section 'Example  %
% raw data file', Supplement).                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%============================================================================================================================================

function [r_Q1, r_Q2, r_Q3, r_Q1_Q3_S_arranged, r_Q3_resampl, annot_T_S_Q2, annot_T_S_Q3, D_str] = ResponseAnalysis(F, annot_W, F_merged, ...
                                                                          annot_W_unique, N_P, annot_P, N_resampl, BarCode, resDir, codeDir)
 
%% Keep Only Response Data (average + replicates)
ind_S_unique = cell2mat(cellfun(@(x)  contains(x, '-S'), annot_W_unique, 'UniformOutput', false));
N_W_S_unique = sum(ind_S_unique);
F_merged_S = F_merged(ind_S_unique, :);
annot_W_S_unique = annot_W_unique(ind_S_unique);

ind_S = cell2mat(cellfun(@(x)  contains(x, '-S'), annot_W, 'UniformOutput', false));
N_W_S = sum(ind_S);
F_S = F(ind_S, :);
annot_W_S = annot_W(ind_S);
                                                                                
%% Intra-Plate Averaging & Random Leave-Ins of Intra-Plate Replicates for Stability Check
if (N_W_S_unique ~= N_W_S)
     
    %%%%%%% Random Leaving-In %%%%%%%
    F_S_resampl = cell(1, N_resampl);
    for j = 1 : N_resampl
        
        F_S_resampl_tmp = zeros(N_W_S_unique, N_P);
        for i = 1 : N_W_S_unique
           
            ind = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W_S, repmat(cellstr(annot_W_S_unique{i}), N_W_S, 1), 'UniformOutput', false));
            F_S_resampl_tmp_ind = F_S(ind, :);
            N_tmp = sum(ind);
            
            if (N_tmp > 1)
                ind_datasample = datasample(1 : N_tmp, N_tmp - 1, 'Replace', false);
                if (N_tmp == 2)
                    F_S_resampl_tmp(i, :) = F_S_resampl_tmp_ind(ind_datasample, :);
                else
                    F_S_resampl_tmp(i, :) = median(F_S_resampl_tmp_ind(ind_datasample, :));
                end
            else
                F_S_resampl_tmp(i, :) = F_S_resampl_tmp_ind;
            end
            
        end
        F_S_resampl{j} = F_S_resampl_tmp;
   
    end

end  

%% Visualization of Response Data after Quality Control and Intra-Plate Averaging
figure();
h = heatmap(annot_P, annot_W_S_unique, round(F_merged_S, 3, 'significant'));
h.Title = 'Pre-Processed Raw Response Data (intra-plate median merging)';
h.XLabel = 'proteins';
h.YLabel = 'wells';
h.Colormap = cool;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat('Pre_Processed_Raw_Response_Data_Median_', BarCode, '.eps'));
clear tmp;
cd(codeDir);

%% Identify string used for the disease asscociated cells
annot_W_S_unique_split = cellfun(@(x) strsplit(x, '-'), annot_W_S_unique, 'UniformOutput', false);
annot_W_S_unique_split = cellfun(@(x) x{1}, annot_W_S_unique_split, 'UniformOutput', false);
state_ID = unique(annot_W_S_unique_split, 'stable');
no_state_ID = length(state_ID);
if (no_state_ID ~= 2)
    
    error('Incompatible cell states.');
    
end
state_ID_H_ind = cell2mat(cellfun(@(x, y) strcmpi(x, y), state_ID, repmat({'H'}, no_state_ID, 1), 'UniformOutput', false));
H_str = state_ID{state_ID_H_ind, 1};
D_str = state_ID{~state_ID_H_ind, 1};

%% Question Q1: normalized release differences between D,To,Sy and H,To,Sy cells
ind_D_UT_S = cell2mat(cellfun(@(x, y) contains(x, y), annot_W_S_unique, repmat(cellstr(strcat(D_str, '-UT-')), N_W_S_unique, 1), ...
                                                                                                                        'UniformOutput', false));
ind_H_UT_S = cell2mat(cellfun(@(x, y) contains(x, y), annot_W_S_unique, repmat(cellstr(strcat(H_str, '-UT-')), N_W_S_unique, 1), ...
                                                                                                                          'UniformOutput', false));

annot_ind_D_UT_S = annot_W_S_unique(ind_D_UT_S);
F_D_UT_S = F_merged_S(ind_D_UT_S, :);
annot_ind_H_UT_S = annot_W_S_unique(ind_H_UT_S);
F_H_UT_S = F_merged_S(ind_H_UT_S, :);

annot_ind_D_UT_S_split = cellfun(@(x) strsplit(x, '-'), annot_ind_D_UT_S, 'UniformOutput', false);

% Unique Stimulations/Provocations
annot_ind_D_UT_S_split = cellfun(@(x) x{end}, annot_ind_D_UT_S_split, 'UniformOutput', false);

% Number of Unique Stimulations/Provocations
N_S_Q = length(annot_ind_D_UT_S_split);

% Rearrange all corresponding elements in the healthy group so that the
% stimulations match
F_H_UT_S_matched = zeros(N_S_Q, N_P);
annot_S_X = cell(N_S_Q, 1);
for i = 1 : N_S_Q

    S_tmp = annot_ind_D_UT_S_split{i};
    ind = cell2mat(cellfun(@(x, y) contains(x, y), annot_ind_H_UT_S, repmat(cellstr(S_tmp), N_S_Q, 1), 'UniformOutput', false));
    if (sum(ind) == 0) % if not a match then no comparison is possible

        F_D_UT_S(i, :) = []; % remove the instance from the disease group

    else

        F_H_UT_S_matched(i, :) = F_H_UT_S(ind, :);
        annot_S_X{i, 1} = S_tmp;

    end

end

ind_isempty = cell2mat(cellfun(@(x) isempty(x), annot_S_X, 'UniformOutput', false));
annot_S_X = annot_S_X(~ind_isempty);

r_Q1 = NormalizeDiff(F_D_UT_S, F_H_UT_S_matched);

figure();
h = heatmap(annot_P, annot_S_X, round(r_Q1, 3, 'significant'));
h.Title = sprintf('Normalized release differences: %s,UT,S vs. H,UT,S (Therapeutic Need)', D_str);
h.XLabel = 'proteins';
h.Colormap = cool;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat('Q1_S_', BarCode, '.eps'));
clear tmp;
cd(codeDir);

%% Question Q2: normalized release differences between D,Tx,Sy and D,To,Sy cells
ind_D_T_S = cell2mat(cellfun(@(x, y) contains(x, y), annot_W_S_unique, repmat(cellstr(strcat(D_str, '-T')), N_W_S_unique, 1), ...
                                                                                                        'UniformOutput', false));

% Retrieve Data for D-T-S (annotations + fluorescence values)
annot_ind_D_T_S = annot_W_S_unique(ind_D_T_S);
F_D_T_S = F_merged_S(ind_D_T_S, :);

% Retrieve Treatments for D-T-S
annot_T_split = cellfun(@(x) strsplit(x, '-'), annot_ind_D_T_S, 'UniformOutput', false);
annot_T_split = cellfun(@(x) x{2}, annot_T_split, 'UniformOutput', false);

% Unique Treatments
annot_T_split_unique = unique(annot_T_split, 'stable');

% Number of Unique Treatments
N_T = length(annot_T_split_unique);

% Retrieve Data for D-UT-S
F_D_UT_S = F_merged_S(ind_D_UT_S, :);

% Number of Unique Stimulations/Provocations
N_S_Q = length(annot_ind_D_UT_S_split);

annot_T_S_Q2 = cell(N_S_Q, 1);
r_Q2 = cell(N_S_Q, 1);
for j = 1 : N_S_Q
    
    r_Q2_S = zeros(N_T, N_P);
    annot_T_S_Q2_S = cell(N_T, 1);
    for i = 1 : N_T
        
        annot_T_S_tmp = strcat(annot_T_split_unique{i}, '-', annot_ind_D_UT_S_split{j});
        ind = cell2mat(cellfun(@(x, y) contains(x, y), annot_ind_D_T_S, repmat(cellstr(annot_T_S_tmp), length(annot_ind_D_T_S), 1), ...
                                                                                                               'UniformOutput', false));
        if (sum(ind) ~= 0)

            r_Q2_tmp = NormalizeDiff(F_D_T_S(ind, :), F_D_UT_S(j, :));
            annot_T_S_Q2_S{i} = annot_T_S_tmp;
            r_Q2_S(i, :) = r_Q2_tmp;

        end
        
    end
    
    ind_isempty = cell2mat(cellfun(@(x) isempty(x), annot_T_S_Q2_S, 'UniformOutput', false));
    annot_T_S_Q2_S = annot_T_S_Q2_S(~ind_isempty);
    r_Q2_S = r_Q2_S(~ind_isempty, :);
    
    annot_T_S_Q2b_to_sort = cellfun(@(x) strsplit(x, '-'), annot_T_S_Q2_S, 'UniformOutput', false);
    annot_T_S_Q2b_to_sort = cellfun(@(x) x{1}, annot_T_S_Q2b_to_sort, 'UniformOutput', false);
    annot_T_S_Q2b_to_sort = cell2mat(cellfun(@(x) str2double(x(2:end)), annot_T_S_Q2b_to_sort, 'UniformOutput', false));
    [~, sort_ind] = sort(annot_T_S_Q2b_to_sort, 'ascend');
    
    figure();
    h = heatmap(annot_P, annot_T_S_Q2_S(sort_ind), round(r_Q2_S(sort_ind, :), 3, 'significant'));
    h.Title = sprintf('Normalized release differences: %s,T,%s vs. %s,UT,%s (Modulation Capacity)', D_str, annot_ind_D_UT_S_split{j}, ...
                                                                                                                D_str, annot_ind_D_UT_S_split{j});
    h.XLabel = 'proteins';
    h.Colormap = cool;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    cd(resDir);
    print -depsc -painters -r400 tmp
    movefile('tmp.eps', strcat('Q2_', annot_ind_D_UT_S_split{j}, '_', BarCode, '.eps'));
    clear tmp;
    cd(codeDir);
    
    annot_T_S_Q2{j} = annot_T_S_Q2_S;
    r_Q2{j} = r_Q2_S;
    
end

%% Question Q3: normalized release differences between D,Tx,Sy and H,To,Sy cells
    
% Find unique stimulations based on H,To,Sy cells
annot_ind_H_UT_S_split = cellfun(@(x) strsplit(x, '-'), annot_ind_H_UT_S, 'UniformOutput', false);
annot_ind_H_UT_S_split = cellfun(@(x) x{end}, annot_ind_H_UT_S_split, 'UniformOutput', false);

% Number of Unique Stimulations/Provocations
N_S_Q = length(annot_ind_H_UT_S_split);

annot_T_S_Q3 = cell(N_S_Q, 1);
r_Q3 = cell(N_S_Q, 1);
for j = 1 : N_S_Q
    
    r_Q3_S = zeros(N_T, N_P);
    annot_T_S_Q3_S = cell(N_T, 1);
    for i = 1 : N_T
        
        annot_T_S_tmp = strcat(annot_T_split_unique{i}, '-', annot_ind_H_UT_S_split{j}); % form string consisting of treatment + stimuli name
        ind = cell2mat(cellfun(@(x, y) contains(x, y), annot_ind_D_T_S, repmat(cellstr(annot_T_S_tmp), length(annot_ind_D_T_S), 1), ...
                                                                                                                        'UniformOutput', false));
        
        if (sum(ind) ~= 0)

            r_Q3_tmp = NormalizeDiff(F_D_T_S(ind, :), F_H_UT_S(j, :));
            annot_T_S_Q3_S{i, 1} = annot_T_S_tmp;
            r_Q3_S(i, :) = r_Q3_tmp;

        end
        
    end
    
    ind_isempty = cell2mat(cellfun(@(x) isempty(x), annot_T_S_Q3_S, 'UniformOutput', false));
    annot_T_S_Q3_S = annot_T_S_Q3_S(~ind_isempty);
    r_Q3_S = r_Q3_S(~ind_isempty, :);
    
    annot_T_S_Q2c_to_sort = cellfun(@(x) strsplit(x, '-'), annot_T_S_Q3_S, 'UniformOutput', false);
    annot_T_S_Q2c_to_sort = cellfun(@(x) x{1}, annot_T_S_Q2c_to_sort, 'UniformOutput', false);
    annot_T_S_Q2c_to_sort = cell2mat(cellfun(@(x) str2double(x(2:end)), annot_T_S_Q2c_to_sort, 'UniformOutput', false));
    [~, sort_ind] = sort(annot_T_S_Q2c_to_sort, 'ascend');
    
    figure();
    h = heatmap(annot_P, annot_T_S_Q3_S(sort_ind), round(r_Q3_S(sort_ind, :), 3, 'significant'));
    h.Title = sprintf('Normalized release differences: %s,T,%s  vs. H,UT,%s (Restoration Capacity)', D_str, annot_ind_H_UT_S_split{j}, ...
                                                                                                                annot_ind_H_UT_S_split{j});
    h.XLabel = 'proteins';
    h.Colormap = cool;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    cd(resDir);
    print -depsc -painters -r400 tmp
    movefile('tmp.eps', strcat('Q3_', annot_ind_H_UT_S_split{j}, '_', num2str(BarCode), '.eps'));
    clear tmp;
    cd(codeDir);
    
    annot_T_S_Q3{j} = annot_T_S_Q3_S;
    r_Q3{j} = r_Q3_S;
    
end

r_Q1_Q3_S_arranged  = zeros(N_S_Q, N_P);
for i = 1 : N_S_Q
   
    ind_o = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_S_X, repmat(cellstr(annot_ind_H_UT_S_split{i}), length(annot_S_X), 1), ...
                                                                                                           'UniformOutput', false));
    r_Q1_Q3_S_arranged(i, :) = r_Q1(ind_o, :);
    
end

%% Ratios r_Q2c for all different stimulations and validation datasets
r_Q3_resampl = cell(N_S_Q, N_resampl);
for k = 1 : N_resampl
   
    F_S_resampl_tmp = F_S_resampl{k};
    F_D_T_S_resampl = F_S_resampl_tmp(ind_D_T_S, :);
    F_H_UT_S_resampl = F_S_resampl_tmp(ind_H_UT_S, :);
    
    for j = 1 : N_S_Q
       
        r_Q3_resampl_tmp = zeros(N_T, N_P);
        annot_T_S_tmp = annot_T_S_Q3{j};
        for i = 1 : N_T
           
            ind = cell2mat(cellfun(@(x, y) contains(x, y), annot_ind_D_T_S, repmat(cellstr(annot_T_S_tmp{i}), length(annot_ind_D_T_S), 1), ...
                                                                                                                      'UniformOutput', false));
            if (sum(ind) ~= 0)

                tmp = NormalizeDiff(F_D_T_S_resampl(ind, :), F_H_UT_S_resampl(j, :));
                r_Q3_resampl_tmp(i, :) = tmp;

            end
            
        end
        r_Q3_resampl{j, k} = r_Q3_resampl_tmp;
        
    end

end

                                                                                                          
end

%============================================================================================================================================

