% DESCRIPTION
%   Divide all V2V links into numC clusters in an attempt to generate
%   minimum intra-cluster interference: MAX K-CUT based on Table III of
%   
%   R. Y. Chang, et al., "multicell OFDMA downlink resource allocation
%   using a grpah framework,"  IEEE Trans. Veh. Technol. vol. 58, no. 7,
%   pp. 3494-3507, Sep. 2009. 
%
% INPUT
%   alpha_kk_: numV2V x numV2V, the interference levels between V2V links
%   numC: the number of clusters
%
% OUTPUT
%   clusterMat: numC x numV2V/numC, clustering result, where each row
%   dentoes a cluster with each element representing the V2V index.
%
% Le Liang, Georgia Tech, Aug. 27, 2017


function [ clusterMat ] = genClusterMaxCut(alpha_kk_, numC)

numV2V = size(alpha_kk_,2);
clusterMat = zeros(numC, numV2V); % store V2V indices (used in indV2Vtx/indV2Vrx) 
                                  % of each cluster (row)

indPerm = randperm(numV2V);
tmp = indPerm(1:numC);

% arbitrarily assign one link to each cluster
for ii = 1:numC
    clusterMat(ii,tmp(ii)) = 1;
end

indV2V = 1:numV2V;
indV2V(indPerm(1:numC)) = [];

for ii = 1:(numV2V-numC)
    % find the smallest delta over all clusters
    delta = zeros(numC,1);
    for ic = 1:numC
        for kk = 1:numV2V
            if clusterMat(ic,kk) == 1
                delta(ic) = delta(ic) + alpha_kk_(kk,indV2V(ii)) + alpha_kk_(indV2V(ii),kk);
            end
        end
    end
    
    [~, ind] = min(delta);
    clusterMat(ind, indV2V(ii)) = 1;
end


end




