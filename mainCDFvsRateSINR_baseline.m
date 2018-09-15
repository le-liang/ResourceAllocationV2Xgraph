% DESCRIPTION:
%   Compare the cummulative distribution function (CDF) of V2I link
%   capacity and V2V SINR when using different resource allocation schemes
%
% Le Liang, Georgia Tech, Feb 22, 2018
% Updated

% Part V: the baseline algorithm

tic

rndseed0 = 2000;
rndseed = 2000;

numCDF = 1e3; % number of fast fading realizations
rng(rndseed0); % control the random seed for randn, randi, rand

%% Parameters setup
infty = 2000; % used as infinity in the simulation
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
v = 70; % velocity
d_avg = 2.5.*v/3.6; % average inter-vehicle distance according to TR 36.885

% QoS parameters for V2I and V2V links
dB_gamma0 = 5; % SINR_min for V2V in dB
p0 = 0.01; % outage probability for V2V
dB_sig2 = -114; % noise power in dBm
dB_gamma_slow = getSlowSINR(dB_gamma0,p0); % SINR in slow CSI

%%%% dB to linear scale conversion
sig2 = 10^(dB_sig2/10);
gamma0 = 10^(dB_gamma0/10);
Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);
gammaSlow = 10^(dB_gamma_slow/10);
gammaProp = gamma0/(-log(1-p0));

numV2I = 10;
numV2V = 3*numV2I;

sumRateBase = zeros(numCDF,1);
sinrBase = zeros(numCDF,1);
sumRateSun = zeros(numCDF,1);
%%
while (1)
    infeasible = 0;
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
    
    %% CROWN scheme in Sun2016Cluster
    numC = 10; % number of clusters, divisible by numV2V
    % divide V2V links into numC clusters, with each row representing a
    % cluster and the elements denoting V2V indices.
    clusterMat = genClusterSun(alpha_kk_,numC);
    x_mk = zeros(numV2I,numV2V); % spectrum sharing indicator
    Pd_mk = zeros(numV2I,numV2V); % V2V TX power
    Pc_m = zeros(numV2I,numV2V); % V2I TX power
    for ic = 1:numC % cycle through each V2V cluster
        % compute rateTmp for numV2I x numV2I combinations in the cluster
        rateTmp = zeros(numV2I,numV2I);
        indCurrentV2V = clusterMat(ic,:); % index of V2V in the current cluster
        optPd = zeros(numV2I,numV2V,numV2I); % optimal powers for V2V
        optPc = zeros(numV2I,numV2V,numV2I);
        for m = 1:numV2I
            for k = 1:numV2I
                if k > length(indCurrentV2V) % dummy vertices
                    alpha_mB = alpha_mB_(m);
                    rateTmp(m,k) = log2(1 + Pc_max*alpha_mB/sig2);
                    optPd(m,:,k) = zeros(1,numV2V);
                    optPc(m,:,k) = Pc_max;
                else
                    indV2Vm = find(x_mk(m,:)==1); % index set of V2Vs already on the V2I m
                    indV2Vm = [indV2Vm, indCurrentV2V(k)];
                    
                    alpha_m = alpha_mk_(m,indV2Vm);
                    alpha_kk = alpha_kk_(indV2Vm,indV2Vm);
                    alpha_mB = alpha_mB_(m);
                    alpha_kB = alpha_kB_(indV2Vm);
                    % Compute optimal V2I and V2Vs power
                    [rate, Pc_opt, Pd_opt] = calOptPower(alpha_m(:), ...
                        alpha_kk, alpha_mB, alpha_kB, gammaSlow, sig2, Pd_max, Pc_max);
                    
                    if rate < 0
                        rateTmp(m,k) = -infty;
                        continue; % continue to the next V2V-V2I pair
                    end
                    
                    optPd(m,indV2Vm,k) = Pd_opt;
                    optPc(m,indV2Vm,k) = Pc_opt;
                    rateTmp(m,k) = rate;
                end
            end
        end
        % Hungarian algorithm
        [assignCurrentCluster, cost] = munkres(-rateTmp);
        if cost > 0
            % the original problem is infeasible
            infeasible = 1;
            break;
        end
        for ii = 1:numV2I
            % exclude dummy vertices
            if assignCurrentCluster(ii) <= length(indCurrentV2V)
                % Refresh the spectrum sharing indicator matrix
                x_mk(ii,indCurrentV2V(assignCurrentCluster(ii))) = 1;
                
                % Refresh the optimal power controls for all V2Vs and V2I
                Pd_mk(ii,:) = optPd(ii,:,assignCurrentCluster(ii));
                Pc_m(ii,:) = optPc(ii,:,assignCurrentCluster(ii));
            end
        end
    end
    
    %%
    if infeasible == 1
        continue;
    end
    % The proposed two algorithms are always better than CROWN, thus this
    % is enough for infeasibility reporting. 
    
    %% Proposed baseline algorithm   
    clusterMatProp = genClusterMaxCut(alpha_kk_,numV2I); %  
    rng(rndseed);
    iChann = 1;
    while iChann <= numCDF
        h_mBf_ = abs(randn(numV2I,numV2I) + 1j*randn(numV2I,numV2I)).^2/2;
        h_kBf_ = abs(randn(numV2V,numV2I) + 1j*randn(numV2V,numV2I)).^2/2;
        % simple algorithm
        [sumV2Isimple, x_mfkSimple, Pd_mfkSimple, Pc_mfkSimple] = solveMatching3D(clusterMatProp, ...
            alpha_mB_, alpha_kB_, alpha_mk_, alpha_kk_, sig2, Pd_max, Pc_max, gammaProp, h_mBf_, h_kBf_);
        
        % CDF of V2I capacity
        sumRateBase(iChann) = sumRateBase(iChann) + sumV2Isimple;
        
        % CDF of V2V reliability
        % sum along f
        x_mk_f_sim = sum(x_mfkSimple,2);
        Pd_mk_f_sim = sum(Pd_mfkSimple,2);
        Pc_mk_f_sim = sum(Pc_mfkSimple,2);
        
        k = 1;
        m = find(x_mk_f_sim(:,k)==1);
        
        Pd = Pd_mk_f_sim(m,x_mk_f_sim(m,:)==1);
        tp = Pc_mk_f_sim(m,x_mk_f_sim(m,:)==1);

        kProSim = find(x_mk_f_sim(m,:)==1);
        g_mk = alpha_mk_(m,kProSim(k)).*abs(randn(1)+1j*randn(1)).^2/2;
        g_kk = alpha_kk_(kProSim,kProSim).*abs(randn(length(kProSim))+1j*randn(length(kProSim))).^2/2;
        
        if isempty(tp)
            continue;
        end
        Pc = tp(1);
        signal = Pd(k)*g_kk(k,k);
        interf = Pc*g_mk + Pd*g_kk(:,k)-signal;
        sinrBase(iChann) = signal/(sig2+interf);
        
        iChann = iChann + 1;
    end
    
    break; % successfully finish one round of slow fading realization
end

%%
figure;
h0 = cdfplot(sumRateBase);
set(h0,'color','k')
set(h0,'linestyle','-')
legend('Baseline')
xlabel('$\sum C_m$ (bps/Hz)', 'interpreter','latex')
ylabel('CDF')
title('')
set(gca,'YScale','log')


figure
sinr = -10:2:50;
marker_0 = zeros(length(sinr),1);
for ind_sinr = 1:length(sinr)
    marker_0(ind_sinr) = sum(sinrBase <= 10^(sinr(ind_sinr)/10))/numCDF;
end
plot(sinr,marker_0,'k--')

grid 
legend('Baseline')
xlabel('SINR (dB)')
ylabel('CDF')
title('')
set(gca,'YScale','log')
ylim([1e-3,1])

%save(sprintf('mainCDFvsRate_5_%d',rndseed))

toc

