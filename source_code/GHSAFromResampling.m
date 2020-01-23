%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  HSAFromResampling function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates and visualizes statistics for the resampling %
% based validation related to the generalized highest single agent(GHSA)%
% analysis. See section 'Resampling statistics' of 'Methods' in the %
% main article text as well as 'GHSA.m' for more details.                %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% HSA_r_Q_resampl: [N_higher_orderxN_resampl] matrix (see below). Each  %
% element (i, j) contains the corresponding GHSA index for combination  %
% treatment i and validation dataset j.                                 %
%                                                                       %
% N_higher_order: number of combination treatments for which GHSA       %
% indices are calculated. The number of rows in the input matrix        %
% HSA_r_Q_resampl is equal to N_higher_order (see above).               %
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
% extra_str_title: extra string to be added in the title of the GHSA    %
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
% Boxplot with the GHSA indices obtained during validation for all      %
% combination treatments.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==============================================================================================================================================

function [] = GHSAFromResampling(HSA_r_Q_resampl, N_higher_order, N_resampl, annot_W_HSA, extra_str_title, filename_str, BarCode, resDir, codeDir)

figure();
boxplot(HSA_r_Q_resampl', 'Notch', 'on', 'Colors', 'm', 'Symbol', 'c+');
grid on;
ylim([-1 1]);
ylabel('GHSA Index', 'FontWeight', 'Bold', 'FontSize', 10, 'FontName', 'Sans Serif');
set(gca, 'Xtick', 1 : N_higher_order, 'XTickLabel', annot_W_HSA, 'FontWeight', 'Bold', 'FontSize', 10, 'FontName', 'Sans Serif');
set(findobj(gca, 'type', 'line'), 'LineWidth', 1.2);
h = findobj(gca, 'Tag','Box');
for j = 1 : length(h)
   patch(get(h(j), 'XData'), get(h(j), 'YData'), 'm', 'FaceAlpha', .1);
end
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'LineWidth', 1.2, 'Color', 'c');


if (isempty(extra_str_title))
    title_str = sprintf('%d resamplings', N_resampl);
else
    title_str = sprintf('%d resamplings, %s', N_resampl, extra_str_title);
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

print -dpdf -painters -r400 -bestfit tmp
movefile('tmp.pdf', strcat(filename_str, '.pdf'));
cd(codeDir);

end

%==============================================================================================================================================

