%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2020  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  extractSubsets function - 20/01/20  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a customized heatmap for a data matrix. It is %
% used in order to produce the heatmaps associated with the raw protein %
% release measurements (see Supplementary Fig. S1-S2).                  %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% M: [Nxd] data matrix, where N is the total number of wells and d the  %
% total number of measured proteins.                                    %
%                                                                       %
% colormap_str: string that defines the colormap to be used. It should  %
% be a valid MATLAB colormap string.                                    %
%                                                                       %
% N_P: number of measured proteins.                                     %
%                                                                       %
% N_W: number of experimental wells.                                    %
%                                                                       %
% annot_P: cell array with as many cells as the number of measured      %
% proteins N_P. A particular cell {i} contains the name of measured     %
% protein i.                                                            %
%                                                                       %
% annot_W: cell array with as many cells as the number of experimental  %
% wells. A particular cell {i} contains the annotation for the cell     %
% state of well i.                                                      %
%                                                                       %
% xtick_fontsize: size for the tick labels on the x-axis.               %
%                                                                       %
% ytick_fontsize: size for the tick labels on the y-axis.               %
%                                                                       %
% title_str: title for the heatmap to be produced.                      %
%                                                                       %
% title_str_fontsize: size for the title of the heatmap.                %
%                                                                       %
% filename_str: filename for saving the heatmap to be produced.         %
%                                                                       %
% resDir: directory where the generated results should be saved.        %
%                                                                       %
% codeDir: directory with the source code.                              %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% Heatmap with all the aforementioned user-defined customizations. NaN  %
% (i.e., missing) values are displayed in white.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==============================================================================================================================================

function [] = generateHeatmap(M, colormap_str, N_P, N_W, annot_P, annot_W, xtick_fontsize, ytick_fontsize, title_str, title_str_fontsize, ...
                                                                                                                filename_str, resDir, codeDir)

figure();
h = imagesc(M);
set(h, 'alphadata', ~isnan(M));
colormap(colormap_str);
colorbar();
set(gca, 'Xtick', 1 : N_P, 'XTickLabel', annot_P, 'FontSize', xtick_fontsize, 'FontWeight', 'Bold', 'FontName', 'Sans Serif');
set(gca, 'YTick', 1 : N_W, 'YTickLabel', annot_W, 'FontWeight', 'Bold', 'FontSize', ytick_fontsize, 'FontName', 'Sans Serif');
xtickangle(90);
title(title_str, 'FontWeight', 'Bold', 'FontSize', title_str_fontsize, 'FontName', 'Sans Serif');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd(resDir);
print -depsc -painters -r400 tmp
movefile('tmp.eps', strcat(filename_str, '.eps'));
clear tmp;
cd(codeDir);


end

%==============================================================================================================================================

