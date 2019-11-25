%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Copyright (C) 2019  Efthymia Chantzi      %%
%%        GNU General Public license v3          %%
%%                 (LICENSE.md)                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  NormalizeDiff function - 19/10/29  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs normalization of protein release differences   %
% per experimental plate, so that values across different plates become %
% comparable. It essentially implements Eq. (1) of the main article,    %
% which is used to answer questions Q1, Q2 and Q3.                      %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% X: raw protein release data stored as a [Nxd] matrix, where N is the  %
% total number of wells corresponding to cell state X and d is the total%
% number of proteins measured. If there is only one well for cell state %
% X, then the input X is a [1xd] row vector.                            %                                          
%                                                                       %
% Y: raw protein release data stored as a [Nxd] matrix, where N is the  %
% total number of wells corresponding to cell state Y and d is the total%
% number of proteins measured. If there is only one well for cell state %
% Y, then the input Y is a [1xd] row vector.                            %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% r: [Nxd] matrix or row vector [1xd] (see above) with the normalized   %
% differences in protein releases as defined by Eq. (1).                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===============================================================================

function r = NormalizeDiff(X, Y)

if (~isvector(X) && isvector(Y))
    
    rows = size(X, 1);
    Y = repmat(Y, rows, 1);
    
elseif (isvector(X) && ~isvector(Y))
    
    rows = size(Y, 1);
    X = repmat(X, rows, 1);
    
end

r = (X - Y)./(X + Y);

end

%===============================================================================


