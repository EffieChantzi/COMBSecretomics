%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  plotHierarchicalKMeans function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates and plots the centroids of the (sub)clusters %
% for the top-down hierarchical K-Means clustering after validation. For%
% more details see section 'Resamplign statistics' of the 'Results' in  %
% the main article.                                                     %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% r_Clusters: cell array {k_optx2} where k_opt is the optimal number of %
% clusters in the 1st hierarchical level. A particular cell array {i, 1}%
% contains any of the ratio matrices r_Q2, r_Q3, as defined by eq. (1), %
% belonging to cluster i in the 1st hierarchical level. A particular    %
% cell array {i, 2} contains a nested {1xK} cell array, where K is the  %
% total number of subclusters for cluster i. A particular nested cell   %
% array {i, 2}{1, k} contains a matrix [N_T_SxN_P], where N_T_S is the  %
% number of treatments in subcluster k of cluster i and N_P denotes the %
% number of measured proteins (see below).                              %
%                                                                       %
% r_Clusters_annot: cell array {k_optx2}, where k_opt is the total      %
% number of different clusters in the 1st hierarchical level. Each cell %
% {i, 1} contains a cell array {Mx1}, where M is the total number of    %
% treatments belonging to cluster i. A particular cell array {i, 2}     %
% contains a nested {1xK} cell array, where K is the total number of    %
% subclusters for cluster i. For each subcluster k the corresponding    %
% treatment annotations are included in a {N_T_Sx1} cell array, where   %
% N_T_S is the number of treatments in subcluster k of cluster i.       % 
%                                                                       %
% r_Q1: row vector containing the ratios defined in question Q1 using   %
% Eq. (1) for all measured proteins. This row vector defines the total  %
% therapeutic need and is used for visualization purposes in order to   %
% ensure that the drug induced protein release differences are not      %
% deviating.                                                            %
%                                                                       % 
% N_P: number of measured proteins.                                     %
%                                                                       %
% annot_P: cell array {1xN_P}. A particular cell {i} contains the anno- %
% tation for the corresponding protein i.                               %
%                                                                       %
% y_label: string for labeling the y-axis.                              %
%                                                                       %
% extra_str_title: extra string to be added in the title of the figure  %
% showing the top-down hierarchical clustering. This is the number of   %
% resamplings to be performed and/or the identifier for the stimulation %
% used (S1, S2, S3 etc). However, this can be modified according to the %
% users' preferences.                                                   %
%                                                                       %
% N_resampl: user-defined number of validation datasets to be created   %
% based on resampling. See section 'Resampling statistics' in the main  %
% article text for more details.                                        %
%                                                                       %
% filename_str: filename for saving the figure to be produced.          %
%                                                                       %
% exhaustive_subset_ID: numeric identifier, which must be set to 1 for  %
% performing exhaustive subset search with the hierarchical clustering  %
% and 0 otherwise.                                                      %
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
% Figure showing the centroids of each (sub)cluster in the form of line %
% graphs. A corresponding legend on the top left corner shows the names %
% of the treatments belonging to the corresponding (sub)cluster with or %
% without the exhaustive subset search. The generated figure is saved in%
% the directory resDir (see above).                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=================================================================================================================================================

function [] = plotHierarchicalKMeans(r_Clusters, r_Clusters_annot, r_Q1, N_P, annot_P, y_label, extra_str_title, N_resampl, filename_str, ...
                                                                                               exhaustive_subset_ID, BarCode, resDir, codeDir)

load 'cool_colormap.mat' 'cool_colormap';
cmap = cool_colormap;
cmap_fig = cell(1, 2);
cmap_fig{1} = cmap(1 : 32, :);
cmap_fig{2} = (cmap(33 : end, :)); 
marker_symbol_array = {'o-', 's-', 'd-', 'h-'};

k_opt = size(r_Clusters, 1);

figure();
hold on;
for r = 1 : k_opt

    count_cmap = 0;
    cmap = cmap_fig{r};   
    marker_symbol = marker_symbol_array{r};
    
    for c = 1 : 2

        r_Clusters_tmp = r_Clusters{r, c};
        r_Clusters_annot_tmp = r_Clusters_annot{r, c};

        if (~isempty(r_Clusters_annot_tmp))

            if (c == 1)
                
                count_cmap = count_cmap + 1;
                cmap_tmp = cmap(count_cmap, :);

                if (size(r_Clusters_tmp, 1) > 1)

                    centroid_val = mean(r_Clusters_tmp);

                else

                    centroid_val = r_Clusters_tmp;

                end

                ind = cell2mat(cellfun(@(x) iscell(x), r_Clusters_annot_tmp, 'UniformOutput', false));
                if (sum(ind) ~= 0)
                    r_Clusters_annot_tmp = [r_Clusters_annot_tmp{:}];
                end
                r_Clusters_annot_tmp_all = concatenateDrugAnnot(r_Clusters_annot_tmp);
               
                % Exhaustive Subset Search Annotation for 1st hierarchical level
                r_Clusters_annot_tmp = cellfun(@(x) x(2 : end), r_Clusters_annot_tmp, 'UniformOutput', false);
                r_Clusters_annot_tmp_all_unique = exhaustiveSubsetSearch(r_Clusters_annot_tmp, '');
                r_Clusters_annot_tmp_all_unique = cellfun(@(x) strcat('T', x), r_Clusters_annot_tmp_all_unique, 'UniformOutput', false);
                r_Clusters_annot_tmp_all_unique = concatenateDrugAnnot(r_Clusters_annot_tmp_all_unique);

                if (exhaustive_subset_ID == 1)
                    plot(1 : N_P, centroid_val, marker_symbol, 'Color', cmap_tmp, 'MarkerFaceColor', cmap_tmp, 'LineWidth', 1.5, 'MarkerSize', 4, ...
                                                                                            'DisplayName', r_Clusters_annot_tmp_all_unique);
                else
                    plot(1 : N_P, centroid_val, marker_symbol, 'Color', cmap_tmp, 'MarkerFaceColor', cmap_tmp, 'LineWidth', 1.5, 'MarkerSize', 4, ...
                                                                                            'DisplayName', r_Clusters_annot_tmp_all);                                                                  
                end


            else

                no_sub = size(r_Clusters_tmp, 2);

                if (no_sub > 0)

                    for i = 1 : no_sub
                        
                        count_cmap = count_cmap + 10;
                        cmap_tmp = cmap(count_cmap, :);

                        centroid_val = r_Clusters_tmp{i};
                        if (~isvector(centroid_val))

                            centroid_val = mean(centroid_val);

                        end

                        r_Clusters_annot_tmp_sub = r_Clusters_annot_tmp{i};
                        ind = cell2mat(cellfun(@(x) iscell(x), r_Clusters_annot_tmp_sub, 'UniformOutput', false));
                        if (sum(ind) ~= 0)
                            r_Clusters_annot_tmp_sub = [r_Clusters_annot_tmp_sub{:}];
                        end
                        r_Clusters_annot_tmp_sub_all = concatenateDrugAnnot(r_Clusters_annot_tmp_sub);
                        
                        % Exhaustive Subset Search Annotation for 2nd hierarchical level
                        r_Clusters_annot_tmp_sub = cellfun(@(x) x(2 : end), r_Clusters_annot_tmp_sub, 'UniformOutput', false);
                        r_Clusters_annot_tmp_sub_all_unique = exhaustiveSubsetSearch(r_Clusters_annot_tmp_sub, '');
                        r_Clusters_annot_tmp_sub_all_unique = cellfun(@(x) strcat('T', x), r_Clusters_annot_tmp_sub_all_unique, ...
                                                                                                                        'UniformOutput', false);
                        r_Clusters_annot_tmp_sub_all_unique = concatenateDrugAnnot(r_Clusters_annot_tmp_sub_all_unique);

                        
                        if (exhaustive_subset_ID == 1)
                            plot(1 : N_P, centroid_val, strcat(marker_symbol, '-'), 'Color', cmap_tmp, 'MarkerFaceColor', cmap_tmp, ...
                                                          'LineWidth', 1.15, 'MarkerSize', 4, 'DisplayName', r_Clusters_annot_tmp_sub_all_unique);
                        else
                            plot(1 : N_P, centroid_val, strcat(marker_symbol, '-'), 'Color', cmap_tmp, 'MarkerFaceColor', cmap_tmp, ...
                                                                'LineWidth', 1.15, 'MarkerSize', 4, 'DisplayName', r_Clusters_annot_tmp_sub_all);
                        end

                    end

                end

            end

        end

    end

end
y_label_split = strsplit(y_label, ' (');
y_label_split = y_label_split{2};
y_label_split = strsplit(y_label_split, ',');
D_str = y_label_split{1};
if (isempty(extra_str_title))
    y_label_r0 = sprintf('%s,UT,US vs. H,UT,US', D_str);
else
    y_label_r0 = sprintf('%s,UT,%s vs. H,UT,%s', D_str, extra_str_title, extra_str_title);
end
plot(1 : N_P, r_Q1, 'p-', 'Color', [.5 .5 .5], 'MarkerFaceColor', [.5 .5 .5], 'LineWidth', 1.0, 'MarkerSize', 4, 'DisplayName', y_label_r0);
hold off;
grid on;
ylim([-1 1]);
xlim([1 N_P]);
ylabel(y_label, 'FontWeight', 'Bold', 'FontName', 'Sans Serif', 'FontSize', 12);
xlabel('proteins', 'FontWeight', 'Bold', 'FontName', 'Sans Serif', 'FontSize', 12);
title_str = sprintf('Hierarchical K-Means Clustering (%d resamplings)', N_resampl);
title(title_str, 'FontWeight', 'Bold', 'FontName', 'Sans Serif', 'FontSize', 12);
l = legend('-DynamicLegend', 'Location', 'NorthWest');
l.FontSize = 10;
l.FontWeight = 'Bold';
l.NumColumns = 2;
set(gca, 'Xtick', 1 : N_P, 'XTickLabel', annot_P, 'FontWeight', 'Bold', 'FontName', 'Sans Serif');
xtickangle(90);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

if (isempty(extra_str_title))
    filename_str = strcat(filename_str, '_', num2str(N_resampl), '_', 'US_', BarCode);
else
    filename_str = strcat(filename_str, '_', num2str(N_resampl), '_', extra_str_title, '_', BarCode);
end

cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat(filename_str, '.eps'));
clear tmp;
cd(codeDir); 
    
end

%=================================================================================================================================================

