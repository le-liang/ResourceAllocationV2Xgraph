
% DESCRIPTION
%   solving the whole 3D matching problem with a given clustering result


function [ sumRate, x_mfk, Pd_mfk, Pc_mfk ] = solveMatching3D(clusterMatProp, alpha_mB_, alpha_kB_, alpha_mk_, alpha_kk_, sig2, Pd_max, Pc_max, gamma, h_mBf_, h_kBf_)
infty = 2e3;
[numV2I,numV2V] = size(alpha_mk_);

rateMFN = zeros(numV2I,numV2I,numV2I);
x_mfk = zeros(numV2I,numV2I,numV2V); % indicator matrix
Pd_mfk = zeros(numV2I,numV2I,numV2V); % opt power for V2V
Pc_mfk = zeros(numV2I,numV2I,numV2V); % opt power for V2I
for f = 1:numV2I
    g_kB_ = alpha_kB_.*h_kBf_(:,f);
    g_mB_ = alpha_mB_.*h_mBf_(:,f);
    for m = 1:numV2I
        for n = 1:numV2I % V2V cluster index
            % indV2V = clusterMatProp(n,:); % indices of V2V in this cluster
            tmp = 1:numV2V;
            indV2V = tmp(clusterMatProp(n,:)>0);
            alpha_m = alpha_mk_(m,indV2V);
            alpha_kk = alpha_kk_(indV2V,indV2V);
            g_mB = g_mB_(m);
            g_kB = g_kB_(indV2V);
            [rate, Pc_opt, Pd_opt] = calOptPower(alpha_m(:), ...
                alpha_kk, g_mB, g_kB, gamma, sig2, Pd_max, Pc_max);
            if rate < 0 % power optimization is infeasible
                rateMFN(m,f,n) = -infty;
                continue;
            end
            rateMFN(m,f,n) = rate;
            
            % store the optimal power values
            Pd_mfk(m,f,clusterMatProp(n,:)>0) = Pd_opt;
            Pc_mfk(m,f,clusterMatProp(n,:)>0) = Pc_opt;
        end
    end
end

small = 0.01;
% rateMFNadj = rateMFN + infty + 1; % avoid negative numbers
rateMFNadj = rateMFN + small;
[assign, ~, ~] = solution3DMatching(rateMFNadj);
if sum(rateMFNadj(assign>0)<0)>0 || length(rateMFNadj(assign>0))<numV2I
% if sum(rateMFNadj(assign>0)<0)>0
    % If the whole layer, e.g., rate(:,:,n)=-infty,
    % solution3DMatching will ignore this layer, making the
    % matching result contain < numV2I elements.
    % If some of rate(m,f,n)=-infty and the matching picks it, the
    % problem is infeasible too.
    sumRate = -1;
    return;
end

sumRate = sum(rateMFN(assign>0));
% refresh x_mfk indicator matrix
for f = 1:numV2I
    for m = 1:numV2I
        for n = 1:numV2I
            if assign(m,f,n)>0
                x_mfk(m,f,clusterMatProp(n,:)>0) = 1;
            end
        end
    end
end
% clear temp powers generated in matching process
Pd_mfk(x_mfk==0) = 0;
Pc_mfk(x_mfk==0) = 0;


