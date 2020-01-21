%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  main script - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main script that the user should run in order to employ   %
% the COMBSecretomics framework.                                        %
%                                                                       %
%                                                                       %
% Upon initiation, a series of user-defined inputs are required; raw    %
% data file, cut-off thresholds (%) for the quality control, the number %
% of resampling based validation datasets to be created and flag for    %
% employing exhaustive subset search during the top-down hierarchical   %
% clustering. Details related to these inputs as well as suggestions on %
% how they should be selected are included in sections 'Example raw data%
% file' and 'User-defined inputs' in the Supplementary Information.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=================================================================================================================================================

clc;
clear;
close all;

fprintf('\nCopyright (C) 2020  Efthymia Chantzi\n GNU GENERAL PUBLIC LICENSE v3');
fprintf('\nComputational Medicine, Dept. of Medical Sciences');
fprintf('\nUppsala University, January 2020\n');

fprintf('\n-------------------------------------------- COMBSecretomics Initiated  --------------------------------------------\n');

codeDir = pwd;

%% Select .csv Specification File(s)
ext = '';
while ((isempty(ext)) || (~strcmp(ext, '.csv')))

    [spec_filename, spec_PathName, ~] = uigetfile('*.csv', 'Select Specification File');
    
    [~, name, ext] = fileparts(spec_filename);
    BarCode = name;

    % Monitor Message for Selected File
    fprintf('\nSpecification File: %s\n', strcat(spec_PathName, spec_filename));
        
    
end
resDir_name = 'Results';
resDir = fullfile(spec_PathName, resDir_name);
if (~exist(resDir, 'dir'))
   
    mkdir(resDir);
    
end

cd(spec_PathName);
T = readtable(spec_filename, 'ReadRowNames', true, 'ReadVariableNames', true);
cd(codeDir);

annot_W = T.Properties.RowNames;        % annot_W: well annotations
annot_W = cellfun(@(x) strsplit(x, '_'), annot_W, 'UniformOutput', false);
annot_W = cellfun(@(x) x{1}, annot_W, 'UniformOutput', false);
N_W = length(annot_W);                  % N_W: number of wells

annot_P = T.Properties.VariableNames;   % annot_P: protein releases
N_P = length(annot_P);                  % N_P: number of protein releases measured

%% Enter cut-off threshold for blank filtering (%)
blank_tau = nan;
while ((isnan(blank_tau)) || (blank_tau <= 0) || ischar(blank_tau))

    blank_tau = str2double(input('\nEnter cut-off threshold for blank filtering (percent):\n', 's'));

end

%% Enter cut-off threshold for coefficient of variation (%)
cv_tau = nan;
while ((isnan(cv_tau)) || (cv_tau <= 0) || ischar(cv_tau))

    cv_tau = str2double(input('\nEnter cut-off threshold for coefficient of variation (percent):\n', 's'));

end

%% Enter number of resampling-based datasets for stability control
N_resampl = nan;
while ((isnan(N_resampl)) || (N_resampl <= 0) || (mod(N_resampl, 1) ~= 0))

    N_resampl = str2double(input('\nEnter number of resamplings:\n', 's'));

end

%% Number of repeats for 2-Level K-Means Clustering to avoid local minima
R = 30;

%% Enter choice for exhaustive subset search
exhaustiveSearch_ID = nan;
while ((isnan(exhaustiveSearch_ID)) || (exhaustiveSearch_ID ~= 1) && (exhaustiveSearch_ID ~= 2))

    exhaustiveSearch_ID = str2double(input('\nPerform exhaustive subset search for K-means clustering?\n 1. Yes\n 2. No\n', 's'));

end

%% Select Analysis
analysis_ID = nan;
while ((isnan(analysis_ID)) || (analysis_ID ~= 1) && (analysis_ID ~= 2))

    analysis_ID = str2double(input('\nSelect analysis mode:\n 1. Unstimulated Cells\n 2. Stimulated Cells\n', 's'));

end

%% Raw Fluorescence Values
F = table2array(T);

%% Quality Control (Pre-Processing)
[F, annot_W, F_merged, annot_W_unique, annot_P, N_P]  = preProcessing(F, N_W, annot_W, N_P, annot_P, blank_tau, cv_tau, BarCode, resDir, codeDir);

%% Main Analyses
time_start = tic;
if (analysis_ID == 1)

    resDir = fullfile(resDir, 'Unstimulated');
    if (~exist(resDir, 'dir'))
        mkdir(resDir);
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Unstimulated Cells %%%%%%%%%%%%%%%%%%%%%%%
    [r_Q1, r_Q2, r_Q3, r_Q3_resampl, N_T, annot_W_US_unique_treatments, D_str] = BaselineAnalysis(F, annot_W, F_merged, annot_W_unique, ...
                                                                                             N_P, annot_P, N_resampl, BarCode, resDir, codeDir);
    
    %% Top-Down Hierarchical K-means              
    ylabel_str = sprintf('Normalized release differences (%s,T,US vs. H,UT,US)', D_str);
    filename_str = 'Q3_K_Means';
 
    % Check Stability of Clustering Using Parts of Data
    r_Clusters_annot_resampl = cell(1, N_resampl);
    fprintf('\n');
    for i = 1 : N_resampl
       
        fprintf('--- Resampling %d of %d ---\n', i, N_resampl);
        r_Clusters_annot_resampl{i} = biHierarchicalClustering(r_Q3_resampl{i}, N_T, annot_W_US_unique_treatments, R, -0.3); 
    end
    fprintf('\n');
    
    r_Clusters_annot_resampl = HierarchicalPartitionResampling(r_Clusters_annot_resampl, N_resampl, '', ...
                                             strcat(filename_str, '_', 'Resampling_', num2str(N_resampl), '_US_', BarCode), resDir, codeDir);
                                                          
    r_Clusters_resampl = KMeansFromResampling(r_Clusters_annot_resampl{1}, r_Q3, N_T, N_P, annot_W_US_unique_treatments);
    
    plotHierarchicalKMeans(r_Clusters_resampl, r_Clusters_annot_resampl{1}, r_Q1, N_P, annot_P, ylabel_str, '', N_resampl, filename_str, ...
                                                                                              exhaustiveSearch_ID, BarCode, resDir, codeDir);
   
    %% Generalized Highest Single Agent Analysis
    [~, annot_W_HSA_US, N_higher_order] = HSA(r_Q3, N_T, annot_W_US_unique_treatments, '', 'Q3', resDir, codeDir, 0);
    filename_str = 'Q3_HSA_Resampling_';
    
    HSA_r_Q3_resampl = zeros(N_higher_order, N_resampl);
    for i = 1 : N_resampl
       
        [HSA_r_Q3_resampl(:, i), ~, ~] = HSA(r_Q3_resampl{i}, N_T, annot_W_US_unique_treatments, '', 'Q3', resDir, codeDir, 0);

    end
    
    HSAFromResampling(HSA_r_Q3_resampl, N_higher_order, N_resampl, annot_W_HSA_US, '', filename_str, BarCode, resDir, codeDir);
   
    %%%%%%%%%%%%%%%%%%%%%% Unstimulated Cells %%%%%%%%%%%%%%%%%%%%%%%
    
else

    resDir = fullfile(resDir, 'Stimulated');
    if (~exist(resDir, 'dir'))
        mkdir(resDir);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%% Stimulated Cells %%%%%%%%%%%%%%%%%%%%%%%
    [r_Q1, r_Q2, r_Q3, r_Q1_Q3_S_arranged, r_Q3_resampl, annot_T_S_Q2, annot_T_S_Q3, D_str] = ResponseAnalysis(F, annot_W, F_merged, ...
                                                                              annot_W_unique, N_P, annot_P, N_resampl, BarCode, resDir, codeDir);
                                                                                                                                                                                               
    %% Top-Down Hierachical K-means 
    %%%%%%% Q3 %%%%%%%
    no_S_Q3 = length(annot_T_S_Q3);
    
    for i = 1 : no_S_Q3 
       
        r_Q3_tmp = r_Q3{i};
        N_T_Q3_tmp = size(r_Q3_tmp, 1);
        annot_T_S_Q3_tmp = annot_T_S_Q3{i};
        annot_T_S_Q3_tmp = cellfun(@(x) strsplit(x, '-'), annot_T_S_Q3_tmp, 'UniformOutput', false);
        S_Q3 = annot_T_S_Q3_tmp{1}{2};
        annot_T_S_Q3_tmp = cellfun(@(x) x{1}, annot_T_S_Q3_tmp, 'UniformOutput', false);
      
        fprintf('\n-------------- Stimulation %s -------------- \n', S_Q3);
        ylabel_str = sprintf('Normalized release differences (%s,T,%s vs. H,UT,%s)', D_str, S_Q3, S_Q3);
        filename_str = 'Q3_K_Means';
        
        
        %% Hierarchical K-means

        % Check Stability of Clustering Using Parts of Data
        r_Clusters_annot_resampl = cell(1, N_resampl);
        fprintf('\n');
        for j = 1 : N_resampl
            fprintf('--- Resampling %d of %d ---\n', j, N_resampl);
            r_Q3_resampl_tmp = r_Q3_resampl{i, j};
            r_Clusters_annot_resampl{j} = biHierarchicalClustering(r_Q3_resampl_tmp, N_T_Q3_tmp, annot_T_S_Q3_tmp, R, -0.3); 
        end
        fprintf('\n');
        
       
        r_Clusters_annot_resampl = HierarchicalPartitionResampling(r_Clusters_annot_resampl, N_resampl, S_Q3, ...
                                       strcat(filename_str, '_', 'Resampling_', num2str(N_resampl), '_', S_Q3, '_', BarCode), resDir, codeDir);
      
        r_Clusters_resampl = KMeansFromResampling(r_Clusters_annot_resampl{1}, r_Q3_tmp, N_T_Q3_tmp, N_P, annot_T_S_Q3_tmp);
        
        plotHierarchicalKMeans(r_Clusters_resampl, r_Clusters_annot_resampl{1}, r_Q1_Q3_S_arranged(i, :), N_P, annot_P, ylabel_str, S_Q3, ...
                                                                    N_resampl, filename_str, exhaustiveSearch_ID, BarCode, resDir, codeDir);   
      
        %% Generalized Highest Single Agent Analysis
        [~, annot_W_HSA_S, N_higher_order] = HSA(r_Q3_tmp, N_T_Q3_tmp, annot_T_S_Q3_tmp, S_Q3, 'Q3', resDir, codeDir, 0);
        filename_str = 'Q3_HSA_Resampling_';
        
        HSA_r_Q3_resampl = zeros(N_higher_order, N_resampl);
        for j = 1 : N_resampl

            r_Q3_resampl_tmp = r_Q3_resampl{i, j};
            [HSA_r_Q3_resampl(:, j), ~, ~] = HSA(r_Q3_resampl_tmp , N_T_Q3_tmp, annot_T_S_Q3_tmp, S_Q3, 'Q3', resDir, codeDir, 0);

        end
        
        HSAFromResampling(HSA_r_Q3_resampl, N_higher_order, N_resampl, annot_W_HSA_S, S_Q3, filename_str, BarCode, resDir, codeDir);
   
    end                                                                                                              
    %%%%%%%%%%%%%%%%%%%%%% Stimulated Cells %%%%%%%%%%%%%%%%%%%%%%%
    
end


time_end = toc(time_start);
t = datevec(time_end./(60*60*24));
fprintf('\n\n\n-------------------------------------------- COMBSecretomics Completed  --------------------------------------------\n\n');
fprintf('\nTime elapsed: %.2f (hrs) %.2f (min) %.2f (sec)\n\n\n', t(4), t(5), t(6));
BeepSound(1);

%=================================================================================================================================================

