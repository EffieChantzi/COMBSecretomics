
%%  preProcessing function - 19/10/29  %%

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%                 Chantzi Effie                 %%
           %%                COMBSecretomics                %%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs quality control for the collected raw data as  %
% thoroughly described in section 'Quality control explained' in the    %
% Supplementary Information.                                            %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% F: numeric collected dataset (i.e., raw protein releases) for the     %
% whole experimental plate. It is stored in the form of a data matrix,  %
% where rows correspond to experimental wells and columns to measured   %
% proteins.                                                             %
%                                                                       %
% N_W: number of experimental wells included in the input matrix F.     %
%                                                                       %
% annot_W: cell array with as many cells as the number of experimental  %
% wells N_W. A particular cell {i} contains the annotation for the cell %
% state of the corresponding well i.                                    %
%                                                                       %
% N_P: number of measured proteins included in the input matrix F.      %
%                                                                       %
% annot_P: cell array with as many cells as the number of measured      %
% proteins N_P. A particular cell {i} contains the name of the measured %
% protein i.                                                            % 
%                                                                       %
% blank_tau: user-defined cut-off threshold in (%) for the blank        %
% filtering step.                                                       %
%                                                                       %
% cv_tau: user-defined cut-off threshold in (%) for the coefficient of  %
% variation among intra-plate replicate measurements.                   %
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
% F: pre-processed data meaning raw protein releases after all quality  %
% control steps.                                                        %
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
% annot_P: cell array with as many cells as the number of measured      %
% proteins N_P after quality control. A particular cell {i} contains    %
% the name of the corresponding measured protein i.                     %                              
%                                                                       %
% N_P: number of measured proteins remained after quality control.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============================================================================================================================================

function [F, annot_W, F_merged, annot_W_unique, annot_P, N_P]  = preProcessing(F, N_W, annot_W, N_P, annot_P, blank_tau, cv_tau, BarCode, ...
                                                                                                                             resDir, codeDir)

%% All raw protein releases that are smaller than or equal to the 95th percentile 
%  of the corresponding blank wells are set to missing (NaN) values
ind = cell2mat(cellfun(@(x) strcmpi(x, 'blank'), annot_W, 'UniformOutput', false));
N_B = sum(ind);
N_W = N_W - N_B;
annot_W = annot_W(~ind, :);
F_B = F(ind, :);
F = F(~ind, :);
F_B_tau = prctile(F_B, 95);
F_B_tau = repmat(F_B_tau, N_W, 1);
ind = (F <= F_B_tau);       
F(ind) = nan;

%% Visualize all raw protein releases apart from blank wells
path_filename_fig = fullfile(resDir, strcat('Raw_Data_Blank_Filtering_', BarCode, '.eps'));
if (~exist(path_filename_fig, 'file'))
    generateHeatmap(F, 'cool', N_P, N_W, annot_P, annot_W, 8, 5, 'Raw Data (Blank Filtering)', 10, ...
                                                                            sprintf('Raw_Data_Blank_Filtering_%s', BarCode), resDir, codeDir);
end

%% Blank filtering across protein releases
blank_tau = blank_tau/100;
ind_F_P_NaN = isnan(F);
F_P_NaN = sum(ind_F_P_NaN);
tau_P_NaN = round(blank_tau*N_W);
ind_tau_P_NaN = (F_P_NaN >= tau_P_NaN);
N_P_excluded = sum(ind_tau_P_NaN);
if (N_P_excluded > 0)
    
    N_P = N_P - N_P_excluded;
    annot_P = annot_P(~ind_tau_P_NaN);
    F = F(:, ~ind_tau_P_NaN);
    
end

%% Blank filtering across treatments (wells)
ind_F_W_NaN = isnan(F);
F_W_NaN = sum(ind_F_W_NaN, 2);
tau_W_NaN = round(blank_tau*N_P);
ind_tau_W_NaN = (F_W_NaN >= tau_W_NaN);
N_W_excluded = sum(ind_tau_W_NaN);
if (N_W_excluded > 0)
    
    N_W = N_W - N_W_excluded;
    F = F(~ind_tau_W_NaN, :);
    annot_W = annot_W(~ind_tau_W_NaN);
    
end

%% Visualization of raw protein releases after blank filtering
path_filename_fig = fullfile(resDir, strcat('Raw_Data_Blank_Filtering_QC_', BarCode, '.eps'));
if (~exist(path_filename_fig, 'file'))
    generateHeatmap(F, 'cool', N_P, N_W, annot_P, annot_W, 8, 6, 'Raw Data (Blank Filtering, QC)', 10, ...
                                                                            sprintf('Raw_Data_Blank_Filtering_QC_%s', BarCode), resDir, codeDir);                                                                                                
end

%% Imputation of missing (NaN) values
ind_NaN = isnan(F);
if(sum(sum(ind_NaN)) > 0)
    
    F = knnimpute(F')';

    %% Visualization of raw protein releases after imputation
    path_filename_fig = fullfile(resDir, strcat('Raw_Data_Imputation_', BarCode, '.eps'));
    if (~exist(path_filename_fig, 'file'))
        generateHeatmap(F, 'cool', N_P, N_W, annot_P, annot_W, 8, 6, 'Raw Data (Imputation)', 10, sprintf('Raw_Data_Imputation_%s', BarCode), ...
                                                                                                                                resDir, codeDir);
    end
                                                           
end

%% Coefficient of variation filtering
annot_W_unique = unique(annot_W, 'stable');
no_unique_annot_W = length(annot_W_unique);
if (no_unique_annot_W ~= N_W)
    
    F_merged = zeros(no_unique_annot_W, N_P);
    F_CV = zeros(no_unique_annot_W, N_P);
    
    for i = 1 : no_unique_annot_W
        
        ind = cell2mat(cellfun(@(x, y) strcmp(x, y), annot_W, repmat(cellstr(annot_W_unique{i}), N_W, 1), 'UniformOutput', false));
        F_ind = F(ind, :);
        no_intra_repli = sum(ind);
        if(no_intra_repli > 1)
            F_merged(i, :) = median(F_ind);
        else
            F_merged(i, :) = F_ind;
        end

        std_ind = std(F_ind);
        mean_ind = mean(F_ind);
        F_CV(i, :) = std_ind./mean_ind;
        
    end
    
end

%% Visualization of coefficient of variation among protein releases
F_CV_median = round(median(F_CV), 3, 'significant');
path_filename_fig = fullfile(resDir, strcat('Raw_Data_CV_Proteins_Median_', BarCode, '.eps'));
if (~exist(path_filename_fig, 'file'))
    figure();
    h = heatmap(annot_P, {'median CV'}, F_CV_median);
    h.Title = 'Intra-Plate Median Coefficient of Variation';
    h.XLabel = 'proteins';
    h.Colormap = cool;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    cd(resDir);
    print -depsc -painters -r400 tmp
    movefile('tmp.eps', strcat('Raw_Data_CV_Proteins_Median_', BarCode, '.eps'));
    clear tmp;
    cd(codeDir);
end

%% Cut-off proteins for which the coefficient of variation is bigger than 
%  the corresponding user-defined threshold
cv_tau = cv_tau/100;
ind_high_CV = (F_CV_median > cv_tau);

F(:, ind_high_CV) = [];
F_merged(:, ind_high_CV) = [];
N_P = N_P - sum(ind_high_CV);
annot_P(ind_high_CV) = [];

%% CSV file with raw protein releases after all quality control steps (i.e., pre-processed data)
path_filename_fig = fullfile(resDir, strcat('Pre_Processed_Raw_Data_', BarCode, '.csv'));
if (~exist(path_filename_fig, 'file'))
    F_table = num2cell(F);
    F_table = [annot_W F_table];
    annot_P_table = [{''} annot_P];
    F_table = [annot_P_table ; F_table];
    cd(resDir);
    writetable(cell2table(F_table), strcat('Pre_Processed_Raw_Data_', BarCode, '.csv'), 'WriteVariableNames', 0);
    cd(codeDir);
end

%% Visualization of raw protein releases after all quality control steps (i.e., pre-processed data) 
% and intra-plate averaging of replicate measurements 
path_filename_fig = fullfile(resDir, strcat('Pre_Processed_Raw_Data_Median_', BarCode, '.eps'));
if (~exist(path_filename_fig, 'file'))
    figure();
    h = heatmap(annot_P, annot_W_unique, round(F_merged, 3, 'significant'));
    h.Title = 'Pre-Processed Raw Data (intra-plate median merging)';
    h.XLabel = 'proteins';
    h.Colormap = cool;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    cd(resDir);
    print -depsc -painters -r400 tmp
    movefile('tmp.eps', strcat('Pre_Processed_Raw_Data_Median_', BarCode, '.eps'));
    clear tmp;
    cd(codeDir);
end

end

%============================================================================================================================================

