%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2019  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  HSAFromResampling function - 19/10/29  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates and visualizes statistics for the resampling %
% based validation related to the generalized highest single agent (HSA)%
% analysis. See section 'Resampling statistics' of the 'Results' in the %
% main article text as well as 'HSA.m' for more details.                %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% HSA_r_Q_resampl: [N_higher_orderxN_resampl] matrix (see below). Each  %
% element (i, j) contains the corresponding HSA index for combination   %
% treatment i and validation dataset j.                                 %
%                                                                       %
% N_higher_order: number of combination treatments for which HSA indices%
% are calculated. The number of rows in the input matrix HSA_r_Q_resampl%
% is equal to N_higher_order (see above).                               %
%                                                                       %
% N_resampl: user-defined number of validation datasets to be created   %
% based on resampling. See section 'Resampling statistics' in the main  %
% article text for more details.                                        %
%                                                                       %
% annot_W_HSA: annotations for combination treatments as a cell array   %
% {N_higher_orderx1} (see above). A particular cell {i} contains the    %
% corresponding annotation for the combination treatment in well i. See %
% function 'HSA' for more details.                                      %
%                                                                       %
% extra_str_title: extra string to be added in the title of the HSA     %
% figure. This is the identifier for the stimulation used (S1, S2, S3   %
% etc) and/or the number of resamplings However, this can be modified   %
% according to the users' preferences.                                  %
%                                                                       %
% filename_str: filename for the figure to be saved.                    %
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
% Heatmap with descriptive statistics (i.e., 2.5^th, 50^th and 97.5^th  %
% percentiles) for the HSA indices obtained during validation for all   %
% combination treatments. The color coding of each row is based on the  %
% corresponding median value (50^th percentile). The generated heatmap  %
% is saved in the directory resDir (see above).                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==============================================================================================================================================

function [] = HSAFromResampling(HSA_r_Q_resampl, N_higher_order, N_resampl, annot_W_HSA, extra_str_title, filename_str, BarCode, resDir, codeDir)

p_low = round(prctile(HSA_r_Q_resampl, 2.5, 2), 2, 'significant');
p_median = round(prctile(HSA_r_Q_resampl, 50, 2), 2, 'significant');
p_high = round(prctile(HSA_r_Q_resampl, 97.5, 2), 2, 'significant');
p_low_strings = num2cell(p_low);
p_low_strings = cellfun(@(x) num2str(x), p_low_strings, 'UniformOutput', false);
p_median_strings = num2cell(p_median);
p_median_strings = cellfun(@(x) num2str(x), p_median_strings, 'UniformOutput', false);
p_high_strings = num2cell(p_high);
p_high_strings = cellfun(@(x) num2str(x), p_high_strings, 'UniformOutput', false);
p_strings = cellfun(@(x, y, z) strcat('[', x, ',', y, ',', z, ']'), p_low_strings, p_median_strings, p_high_strings, 'UniformOutput', false);


figure();
imagesc(p_median);
colormap('cool');
colorbar();
text(ones(N_higher_order, 1), 1 : N_higher_order, p_strings, 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'Bold', ...
                                                                                                                'FontName', 'SansSerif');
set(gca, 'Xtick', 1 : 1, 'XTickLabel', {'[2.5^{th}, 50^{th}, 97.5^{th}] percentiles'}, 'FontSize', 10, 'FontWeight', 'Bold', ...
                                                                                                               'FontName', 'Sans Serif');
set(gca, 'YTick', 1 : N_higher_order, 'YTickLabel', annot_W_HSA, 'FontWeight', 'Bold', 'FontSize', 10, 'FontName', 'Sans Serif');
if (isempty(extra_str_title))
    title_str = sprintf('HSA Index (%d resamplings)', N_resampl);
else
    title_str = sprintf('HSA Index (%d resamplings, %s)', N_resampl, extra_str_title);
end
title(title_str, 'FontWeight', 'Bold', 'FontSize', 10, 'FontName', 'Sans Serif');

if (isempty(extra_str_title))
    filename_str = strcat(filename_str, num2str(N_resampl), '_', 'US_', BarCode);
else
    filename_str = strcat(filename_str, num2str(N_resampl), '_', extra_str_title, '_', BarCode);
end
cd(resDir);
print -depsc -painters -r400 tmp

movefile('tmp.eps', strcat(filename_str, '.eps'));
clear tmp;
cd(codeDir);

end

%==============================================================================================================================================

