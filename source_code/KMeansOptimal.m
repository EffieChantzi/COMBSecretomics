
%%  KMeansOptimal function - 19/10/29  %%

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%                 Chantzi Effie                 %%
           %%                COMBSecretomics                %%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs K-Means clustering by employing MATLAB's built-%
% in function 'kmeans' with euclidean distance as similarity metric. The%
% number clusters to be used is optimized as described in section       %
% 'Mining of extracted response patterns', Materials and Methods of the %
% article:                                                              %
% https://doi.org/10.1186/s12859-019-2908-0 .                           %
%                                                                       %
%                                                                       %
% %%%% INPUTS %%%%                                                      %
% X: input [Nxd] matrix, where N is the number of treatments and d is   %
% the number of protein releases measured per treatment. In terms of the%
% COMBSecretomics framework, X can be the r_Q matrix calculated in Q2 or%
% Q3 when Eq. (1) is employed for the all measured proteins.            %
%                                                                       %
% Kmax: number of maximum number of clusters to be used.                %
%                                                                       %
% repeats: number of repetitions/restarts for K-Means clustering.       %
%                                                                       %
% tau_SSE_drop: drop in the sum of squared errors (SSE) when transitio- %
% ning from k-1 to k as the number of clusters to be used. It must be a %
% negative number (e.g., -0.3, -0.2, -0.5 etc).                         %
%                                                                       %
%                                                                       %
% %%%% OUTPUTS: %%%%                                                    %
% SSE: sum of squared errors as a [repeatsxKmax] matrix (see above). A  %
% particular element (i, j) contains the SSE for repetition i when using%
% j clusters.                                                           %
%                                                                       %
% SSE_drop: drop in SSE (see above) as a [1xKmax-1] matrix. A particular%
% element (i) contains the drop in SSE when transitioning from i to i+1 %
% clusters. For this calculation, the minimum SSE among all repeats for %
% a particular number of clusters is chosen (see above).                %
%                                                                       %
% idx: cell array {repeatsxKmax} (see above). A particular cell {i, j}  %
% contains a column vector with all corresponding cluster indices for   %
% repetition i when using j clusters.                                   %
%                                                                       %
% idx_opt: column vector [Nx1], where N denotes the total number of     %
% treatments to be clustered. It contains the cluster indices for the   %
% optimal number of clusters to be used (i.e., the smallest K that      %
% results in SSE drop more than tau_SSE_drop compared to K-1).          %
%                                                                       %
% k_opt: optimal number of clusters to be used (see above).             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=================================================================================================

function [SSE, SSE_drop, idx, idx_opt, k_opt] = KMeansOptimal(X, Kmax, repeats, tau_SSE_drop)


SSE = zeros(repeats, Kmax);
idx = cell(repeats, Kmax);

for r = 1 : repeats
    
    for k = 1 : Kmax

        SSE_k = 0;
        [idx_tmp, C] = kmeans(X, k);
        idx{r, k} = idx_tmp;
        
        C_vec = unique(idx_tmp, 'stable');  
        numOfC = length(C_vec);
        for i = 1 : numOfC
            
            N_C = sum(idx_tmp == C_vec(i));
            X_C = X(idx_tmp == C_vec(i), :);
            X_Rec = repmat(C(C_vec(i), :), N_C, 1);
            SSE_C = (sum(sum((X_Rec - X_C).^2)));
            SSE_k = SSE_k + SSE_C;

        end
        SSE(r, k) = SSE_k;

    end

end

SSE_min = min(SSE);
SSE_drop = zeros(1, Kmax - 1);
for i = 2 : Kmax
   
    SSE_drop(1, i - 1) = (SSE_min(i) - SSE_min(i - 1))/SSE_min(i - 1);
    
    
end


k_opt_ind = find(SSE_drop <= tau_SSE_drop);

if (~isempty(k_opt_ind))
    
    k_opt = k_opt_ind(1) + 1;
    
else
    
    k_opt_ind = find(SSE_drop == min(SSE_drop));
    if (length(k_opt_ind) == 1)
        k_opt = k_opt_ind(1) + 1;
    else
        k_opt = k_opt_ind + 1;
    end
    
end

SSE_k_opt = SSE(:, k_opt);
idx_opt = idx(:, k_opt);
idx_min_SSE_k_opt = (SSE_k_opt == min(SSE_k_opt));
idx_opt = idx_opt(idx_min_SSE_k_opt);
if (sum(idx_min_SSE_k_opt) > 1)
   
    idx_opt = idx_opt{1};
    
else
    
    if (iscell(idx_opt))
        idx_opt = cell2mat(idx_opt);
    end
    
end


end

%=========================================================================================================

