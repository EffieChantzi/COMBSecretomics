
%%  BaselineAnalysis function - 19/10/29 %%

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%                 Chantzi Effie                 %%
           %%                COMBSecretomics                %%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function normalizes the protein release differences per plate for%
% unstimulated cells. More precisely, this function gives answers to    %
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
% r_Q3_resampl: cell array {1xN_resampl}, where N_resampl is the user-  %
% defined number of validation datasets to be created. A particular cell%
% {i} contains the matrix r_Q3 (see output above) for the corresponding %
% validation dataset i.                                                 %
%                                                                       %
% N_T: number of unique treatments (after quality control and intra-    %
% plate averaging).                                                     %
%                                                                       %
% annot_W_US_unique_treatments: cell array {1xN_T}. A particular cell   %
% {i} contains the annotation of the corresponding treatment i.         %
%                                                                       %
% D_str: string identifier for the disease model of interest. This is   %
% automatically retrieved from the raw data file (see section 'Example  %
% raw data file', Supplement).                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%============================================================================================================================================

function [r_Q1, r_Q2, r_Q3, r_Q3_resampl, N_T, annot_W_US_unique_treatments, D_str] = BaselineAnalysis(F, annot_W, F_merged, annot_W_unique, ...
                                                                                               N_P, annot_P, N_resampl, BarCode, resDir, codeDir)

%% Keep Only Baseline Data (average + replicates)
ind_US_unique = cell2mat(cellfun(@(x)  contains(x, '-US'), annot_W_unique, 'UniformOutput', false));
N_W_US_unique = sum(ind_US_unique);
F_merged_US = F_merged(ind_US_unique, :);
annot_W_US_unique = annot_W_unique(ind_US_unique);

ind_US = cell2mat(cellfun(@(x)  contains(x, '-US'), annot_W, 'UniformOutput', false));
N_W_US = sum(ind_US);
F_US = F(ind_US, :);
annot_W_US = annot_W(ind_US);

%% Intra-Plate Averaging & Random Leave-Ins of Intra-Plate Replicates for Stability Check
if (N_W_US_unique ~= N_W_US)
    
    %%%%%%% Random Leaving-In %%%%%%%
    F_US_resampl = cell(1, N_resampl);
    for j = 1 : N_resampl
        
        F_US_resampl_tmp = zeros(N_W_US_unique, N_P);
        for i = 1 : N_W_US_unique
           
            ind = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W_US, repmat(cellstr(annot_W_US_unique{i}), N_W_US, 1), 'UniformOutput', false));
            F_US_resampl_tmp_ind = F_US(ind, :);
            N_tmp = sum(ind);
            if (N_tmp > 1)
                ind_datasample = datasample(1 : N_tmp, N_tmp - 1, 'Replace', false);
                if (N_tmp == 2)
                    F_US_resampl_tmp(i, :) = F_US_resampl_tmp_ind(ind_datasample, :);
                else
                    F_US_resampl_tmp(i, :) = median(F_US_resampl_tmp_ind(ind_datasample, :));
                end  
            else
                F_US_resampl_tmp(i, :) = F_US_resampl_tmp_ind;
            end
             
        end
        F_US_resampl{j} = F_US_resampl_tmp;
   
    end
   
end

%% Visualization of Baseline Data after QC, Imputation and Intra-Plate Averaging
figure();
h = heatmap(annot_P, annot_W_US_unique, round(F_merged_US, 3, 'significant'));
h.Title = 'Pre-Processed Raw Baseline Data (intra-plate median merging)';
h.XLabel = 'proteins';
h.YLabel = 'wells';
h.Colormap = cool;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat('Pre_Processed_Raw_Baseline_Data_Median_', BarCode, '.eps'));
clear tmp;
cd(codeDir);

%% Identify string used for the disease asscociated cells
annot_W_US_unique_split = cellfun(@(x) strsplit(x, '-'), annot_W_US_unique, 'UniformOutput', false);
annot_W_US_unique_split = cellfun(@(x) x{1}, annot_W_US_unique_split, 'UniformOutput', false);
state_ID = unique(annot_W_US_unique_split, 'stable');
no_state_ID = length(state_ID);
if (no_state_ID ~= 2)
    
    error('Incompatible cell states.');
    
end
state_ID_H_ind = cell2mat(cellfun(@(x, y) strcmpi(x, y), state_ID, repmat({'H'}, no_state_ID, 1), 'UniformOutput', false));
H_str = state_ID{state_ID_H_ind, 1};
D_str = state_ID{~state_ID_H_ind, 1};


%% Question Q1: normalized release differences between D,UT,US and H,UT,US cells
ind_D_UT_US = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W_US_unique, repmat(cellstr(strcat(D_str, '-UT-US')), N_W_US_unique, 1), ...
                                                                                                                        'UniformOutput', false));
ind_H_UT_US = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W_US_unique, repmat(cellstr(strcat(H_str, '-UT-US')), N_W_US_unique, 1), ...
                                                                                                                        'UniformOutput', false));
r_Q1 = NormalizeDiff(F_merged_US(ind_D_UT_US, :), F_merged_US(ind_H_UT_US, :));
r_Q1 = round(r_Q1, 3);

figure();
h = heatmap(annot_P, {'r'}, round(r_Q1, 3, 'significant'));
h.Title = sprintf('Normalized release differences: %s,UT,US vs. H,UT,US (Therapeutic Need)', D_str);
h.XLabel = 'proteins';
h.Colormap = cool;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat('Q1_US_', BarCode, '.eps'));
clear tmp;
cd(codeDir);

%% Question Q2: normalized release differences between D,T,US and D,UT,US cells
ind_D_T_US = cell2mat(cellfun(@(x, y) contains(x, y), annot_W_US_unique, repmat(cellstr(strcat(D_str, '-T')), N_W_US_unique, 1), ...
                                                                                                                    'UniformOutput', false));
N_T = sum(ind_D_T_US);     % number of treatments

r_Q2 = NormalizeDiff(F_merged_US(ind_D_T_US, :), F_merged_US(ind_D_UT_US, :));

annot_W_US_unique_treatments = cellfun(@(x) strsplit(x, '-'), annot_W_US_unique(ind_D_T_US), 'UniformOutput', false);
annot_W_US_unique_treatments = cellfun(@(x) x{2}, annot_W_US_unique_treatments, 'UniformOutput', false);
annot_W_US_unique_treatments_to_sort = cell2mat(cellfun(@(x) str2double(x(2:end)), annot_W_US_unique_treatments, 'UniformOutput', false));
[~, sort_ind] = sort(annot_W_US_unique_treatments_to_sort, 'ascend');
annot_W_US_unique_treatments_sort = annot_W_US_unique_treatments(sort_ind);

figure();
h = heatmap(annot_P, annot_W_US_unique_treatments_sort, round(r_Q2(sort_ind, :), 3, 'significant'));
h.Title = sprintf('Normalized release differences: %s,T,US vs. %s,UT,US (Modulation Capacity)', D_str, D_str);
h.XLabel = 'proteins';
h.Colormap = cool;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat('Q2_US_', BarCode, '.eps'));
clear tmp;
cd(codeDir);

%% Question Q3: normalized release differences between D,T,US and H,UT,US cells
r_Q3 = NormalizeDiff(F_merged_US(ind_D_T_US, :), F_merged_US(ind_H_UT_US, :));

figure();
h = heatmap(annot_P, annot_W_US_unique_treatments_sort, round(r_Q3(sort_ind, :), 3, 'significant'));
h.Title = sprintf('Normalized release differences: %s,T,US  vs. H,UT,US (Restoration Capacity)', D_str);
h.XLabel = 'proteins';
h.Colormap = cool;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat('Q3_US_', BarCode, '.eps'));
clear tmp;
cd(codeDir);

%% Ratios r_Q3 for all different validation datasets
r_Q3_resampl = cell(1, N_resampl);
for i = 1 : N_resampl
   
    F_US_resampl_tmp = F_US_resampl{i};
    r_Q3_resampl{i} = NormalizeDiff(F_US_resampl_tmp(ind_D_T_US, :), F_US_resampl_tmp(ind_H_UT_US, :));
    
end


end

%============================================================================================================================================

