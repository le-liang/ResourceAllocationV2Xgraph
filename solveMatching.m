
% DESCRIPTION
%   solving the whole matching problem with a given clustering result


function [ sumRate, y_mk, Pds_mk, Pcs_m ] = solveMatching(clusterMat, alpha_mB_, alpha_kB_, alpha_mk_, alpha_kk_, sig2, Pd_max, Pc_max, gamma)

infty = 2e3;
[numV2I,numV2V] = size(alpha_mk_);

y_mk = zeros(numV2I, numV2V); % spectrum sharing indicator
Pds_mk = zeros(numV2I,numV2V); % V2V TX power
Pcs_m = zeros(numV2I,numV2V); % V2I TX power
rateMN = zeros(numV2I,numV2I);
for m = 1:numV2I
    for n = 1:numV2I % V2V cluster index
        if sum(clusterMat(n,:)>0)<1
            % empty cluster
            rateMN(m,n) = log2(1 + Pc_max*alpha_mB_(m)/sig2);          
            continue;
        end          
        tmp = 1:numV2V;
        indV2V = tmp(clusterMat(n,:)>0);
        alpha_m = alpha_mk_(m,indV2V);
        alpha_kk = alpha_kk_(indV2V,indV2V);
        g_mB = alpha_mB_(m);
        g_kB = alpha_kB_(indV2V);
        [rate, Pc_opt, Pd_opt] = calOptPower(alpha_m(:), ...
            alpha_kk, g_mB, g_kB, gamma, sig2, Pd_max, Pc_max);
        if rate < 0 % power optimization problem is infeasible
            rateMN(m,n) = -infty;
            continue;
        end
        rateMN(m,n) = rate;
        
        Pds_mk(m,clusterMat(n,:)>0) = Pd_opt;
        Pcs_m(m,clusterMat(n,:)>0) = Pc_opt;
    end
end

% Hungarian algorithm
[assignSlow, cost] = munkres(-rateMN);
if cost > 0
    % the original problem is infeasible
    sumRate = -1;
    return;
end

sumRate = - cost;

for ii = 1:numV2I
    % Refresh the spectrum sharing indicator matrix
    y_mk(ii,clusterMat(assignSlow(ii),:)>0) = 1;
end
% Refresh the optimal power controls for all V2Vs and V2I
Pds_mk(y_mk==0) = 0;
Pcs_m(y_mk==0) = 0;




