% DESCRIPTION
%   Divide all V2V links into numC clusters in an attempt to generate
%   maximum intra-cluster interference according to Algorithm 1 in
%   Sun2016Cluster TWC paper.
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


function [ clusterMat ] = genClusterSun(alpha_kk_, numC)

numV2V = size(alpha_kk_,2);

clusterMat = []; % store V2V indices (used in indV2Vtx/indV2Vrx) of each cluster (row)
residualV2V = 1:numV2V;
for ic = 1:numC
    clusterTmp = [];
    maxTmp = 0;
    k1 = 1;
    k2 = 1;
    for k = 1:length(residualV2V)
        for kk = 1:length(residualV2V)
            if kk ~= k
                if alpha_kk_(residualV2V(k),residualV2V(kk)) >= maxTmp
                    maxTmp = alpha_kk_(residualV2V(k),residualV2V(kk));
                    k1 = k;
                    k2 = kk;
                end
            end
        end
    end
    clusterTmp = [clusterTmp, residualV2V(k1),residualV2V(k2)];
    residualV2V([k1,k2]) = [];
    while length(clusterTmp) < numV2V/numC
        tmpSum = 0;
        % cycle through all elements in residualV2V to find the worst
        % interfering V2V link and put into clusterTmp
        for k = 1:length(residualV2V)
            sumInt = 0;
            for kk = 1:length(clusterTmp)
                sumInt = sumInt + alpha_kk_(residualV2V(k),clusterTmp(kk))...
                    + alpha_kk_(clusterTmp(kk),residualV2V(k));
            end
            if sumInt >= tmpSum;
                V2VtoAdd = k;
                tmpSum = sumInt;
            end
        end
        clusterTmp = [clusterTmp, residualV2V(V2VtoAdd)];
        residualV2V(V2VtoAdd) = [];
    end % end of a cluster
    clusterMat = [clusterMat; clusterTmp];
end
end
