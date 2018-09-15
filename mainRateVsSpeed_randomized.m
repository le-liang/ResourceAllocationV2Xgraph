% DESCRIPTION
%   Plot the sum V2I capacity vs increasing vehicle speed

% By Le Liang, Georgia Tech, Sep. 19, 2017


% function [] = function_mainRateVsSpeed(rndseed0,num,nIter)
tic
num = 1e3;
rndseed0 = 2000;
nIter = 10;

iter = nIter; % number of iterations for greedy/randomized algorithms
channNum = num; % number of fast fading realizations
rng(rndseed0); % control the random seed for randn, randi, rand

%% Parameters setup
dB_Pd_max = 23; % max V2V transmit power in dBm
dB_Pc_max = 23; % max V2I transmit power in dBm

% large scale fading parameters
stdV2V = 3; % shadowing std deviation
stdV2I = 8;

% cell parameter setup
freq = 2; % carrier frequency 2 GHz
radius = 500; % cell radius in meters
bsHgt = 25; % BS height in meters
disBstoHwy = 35; % BS-highway distance in meters
bsAntGain = 8; % BS antenna gain 8 dBi
bsNoiseFigure = 5; % BS noise figure 5 dB

vehHgt = 1.5; % vehicle antenna height, in meters
vehAntGain = 3; % vehicle antenna gain 3 dBi
vehNoiseFigure = 9; % vehicle noise figure 9 dB

numLane = 6;
laneWidth = 4;
v = 20:20:160; % velocity
d_avg_ = 2.5.*v/3.6; % average inter-vehicle distance according to TR 36.885

% QoS parameters for V2I and V2V links
dB_gamma0 = 5; % SINR_min for V2V in dB
p0 = 0.01; % outage probability for V2V
dB_sig2 = -114; % noise power in dBm

%%%% dB to linear scale conversion
sig2 = 10^(dB_sig2/10);
gamma0 = 10^(dB_gamma0/10);
Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);
gammaProp = gamma0/(-log(1-p0));

numV2I = 10;
numV2V = 3*numV2I;
increa = 4; % SA parameter
amplify0 = 1;

%%
sumRateProp = zeros(length(d_avg_),1);
sumRatePropSlow = zeros(length(d_avg_),1);

%%
parfor ind = 1:length(d_avg_) % for parallel computing
    d_avg = d_avg_(ind);
    ich = 1; % channel realization counter
    while ich <= channNum
        %% Generate traffic on the highway
        d0 = sqrt(radius^2-disBstoHwy^2);
        [ Flag, vehPos, indV2I, indV2Vtx, indV2Vrx ] = genVehLinks(d0, laneWidth, numLane, disBstoHwy, d_avg, numV2I, numV2V);
        if Flag == 1
            continue; % generated vehicles are not enough to perform simulation, jump to the next iteration.
        end
        %% random large-scale fading generation
        alpha_mB_ = zeros(numV2I,1);
        alpha_kB_ = zeros(numV2V,1);
        alpha_mk_ = zeros(numV2I,numV2V);
        alpha_kk_ = zeros(numV2V,numV2V);
        for m = 1:numV2I
            % V2I signal link
            dist_mB = norm(vehPos(indV2I(m),:));
            dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, freq) + vehAntGain + bsAntGain - bsNoiseFigure;
            alpha_mB_(m) = 10^(dB_alpha_mB/10);
            
            % Interference from V2I to V2V
            for k = 1:numV2V
                dist_mk = norm(vehPos(indV2I(m),:) - vehPos(indV2Vrx(k),:));
                dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, freq) + 2*vehAntGain - vehNoiseFigure;
                alpha_mk_(m,k) = 10^(dB_alpha_mk/10);
            end
        end
        for k = 1:numV2V
            % Interference from V2V to V2I
            dist_kB = norm(vehPos(indV2Vtx(k),:));
            dB_alpha_kB = genPL('V2I', stdV2I, dist_kB, vehHgt, bsHgt, freq)+ vehAntGain+bsAntGain - bsNoiseFigure;
            alpha_kB_(k) = 10^(dB_alpha_kB/10);
            
            % V2V signal and peer V2V interference
            for kk = 1:numV2V
                dist_kk = norm(vehPos(indV2Vtx(k),:) - vehPos(indV2Vrx(kk),:));
                dB_alpha_kk = genPL('V2V', stdV2V, dist_kk, vehHgt, vehHgt, freq) + 2*vehAntGain - vehNoiseFigure;
                alpha_kk_(k,kk) = 10^(dB_alpha_kk/10);
            end
        end
        
        %% Proposed Slow-CSI case
        clusterMatProp = genClusterMaxCut(alpha_kk_,numV2I); %
        amplify = amplify0;
        for j=1:iter
            for k=1:numV2V
                criteria = zeros(numV2I,1);
                % if only one V2V in the cluster, do not move it out
                if sum(clusterMatProp(clusterMatProp(:,k)==1,:)) == 1
                    continue;
                end
                clusterMatProp(:,k)=0;
                % compute capacity for when putting k into different clusters
                for cl = 1:numV2I
                    clusterMatTmp = clusterMatProp;
                    clusterMatTmp(cl,k) = 1;
                    [ criteria(cl), ~, ~, ~ ] = solveMatching(clusterMatTmp, ...
                        alpha_mB_, alpha_kB_, alpha_mk_, alpha_kk_, sig2, Pd_max, Pc_max, gammaProp);
                end
                
                % randomized algorithm based on criteria
                minRate = min(criteria(criteria>0)); % for numerical stability
                criteriaNew = zeros(numV2I,1);
                criteriaNew(criteria>0) = exp((criteria(criteria>0)-minRate)*amplify);
                criteriaNew = criteriaNew/sum(criteriaNew);
                prob = zeros(1,numV2I);
                for g = 1:numV2I
                    if criteria(g)==-1
                        continue;
                    end
                    prob(g)=sum(criteriaNew(1:g));
                end
                ru = rand(1);
                for g = 1:numV2I
                    if criteria(g)==-1
                        continue;
                    end
                    if ru < prob(g)
                        clusterMatProp(g,k)=1;
                        break;
                    end
                end
            end
            if amplify < 10 % simulated annealing idea
                amplify = amplify*increa;
            end
        end % end of randomized iteration
        
        % Find the optimal powers and spectrum sharing indicator
        [ rateTmp,y_mk,Pds_mk,Pcs_m ] = solveMatching(clusterMatProp, ...
            alpha_mB_, alpha_kB_, alpha_mk_, alpha_kk_, sig2, Pd_max, Pc_max, gammaProp);
        
        if rateTmp < 0
            continue;
        end
        
        %% Proposed algorithm
        clusterMatProp = genClusterMaxCut(alpha_kk_,numV2I); %
        h_mBf_ = abs(randn(numV2I,numV2I) + 1j*randn(numV2I,numV2I)).^2/2;
        h_kBf_ = abs(randn(numV2V,numV2I) + 1j*randn(numV2V,numV2I)).^2/2;
        amplify = amplify0;
        for j=1:iter
            for k=1:numV2V
                criteria = zeros(numV2I,1);
                % if only one V2V in the cluster, do not move it out
                if sum(clusterMatProp(clusterMatProp(:,k)==1,:)) == 1
                    continue;
                end
                clusterMatProp(:,k)=0;
                % compute capacity for when putting k into different clusters
                for cl = 1:numV2I
                    clusterMatTmp = clusterMatProp;
                    clusterMatTmp(cl,k) = 1;
                    [ criteria(cl), ~, ~, ~ ] = solveMatching3D(clusterMatTmp, ...
                        alpha_mB_, alpha_kB_, alpha_mk_, alpha_kk_, sig2, Pd_max, Pc_max, gammaProp,h_mBf_,h_kBf_);
                end
                % randomized algorithm based on criteria
                minRate = min(criteria(criteria>0)); % for numerical stability
                criteriaNew = zeros(numV2I,1);
                criteriaNew(criteria>0) = exp((criteria(criteria>0)-minRate)*amplify);
                criteriaNew = criteriaNew/sum(criteriaNew);
                prob = zeros(1,numV2I);
                for g = 1:numV2I
                    if criteria(g)==-1
                        continue;
                    end
                    prob(g)=sum(criteriaNew(1:g));
                end
                ru = rand(1);
                for g = 1:numV2I
                    if criteria(g)==-1
                        continue;
                    end
                    if ru < prob(g)
                        clusterMatProp(g,k)=1;
                        break;
                    end
                end
            end
            if amplify < 10
                amplify = amplify*increa;
            end
        end
        [sumV2I, ~, ~, ~] = solveMatching3D(clusterMatProp, ...
            alpha_mB_, alpha_kB_, alpha_mk_, alpha_kk_, sig2, Pd_max, Pc_max, gammaProp, h_mBf_, h_kBf_);
        if sumV2I < 0
            continue;
        end
        
        % CDF of V2I capacity
        sumRateProp(ind) = sumRateProp(ind) + sumV2I;
        
        %% generate fast fading for proposed-slow
        for m = 1:numV2I
            kSlow = find(y_mk(m,:)==1);
            g_mB = alpha_mB_(m)*abs((randn(1) + 1j*randn(1))/sqrt(2))^2;
            g_kBslow = alpha_kB_(kSlow).*abs((randn(length(kSlow),1) + 1j*randn(length(kSlow),1))/sqrt(2)).^2;
            
            % proposed slow
            Pc = Pcs_m(m,kSlow);
            signal = Pc(1)*g_mB;
            interf = Pds_mk(m,kSlow)*g_kBslow;
            tmpSlow = log2(1+ signal/(interf + sig2));
            sumRatePropSlow(ind) = sumRatePropSlow(ind) + tmpSlow;
        end
        ich = ich + 1;
        display(ich)
    end   
end

sumRatePropSlow = sumRatePropSlow/channNum;
sumRateProp = sumRateProp/channNum;

figure
plot(v,sumRateProp,'k-+','linewidth',1)
hold on
plot(v, sumRatePropSlow, 'b-*', 'linewidth',1);
grid on
xlabel(Speed, '$v$ (km/h)', 'interpreter','latex')
ylabel('Sum V2I capacity, $\sum\limits_m r_m$ (bps/Hz)', 'interpreter','latex')
legend('P_{max}^c = 23 dBm, Algorithm 6', 'P_{max}^c = 23 dBm, Algorithm 7')

% save(sprintf('mainRateVsSpeed_%d',rndseed0))

toc




