% DESECRIPTION
%   Frist, generate traffic on the highway according to Spatial Poisson Process with
%   density determined by vehicle speed v (represented by avergage
%   inter-vehicle distance d_avg = 2.5*v). 
%   Then, generate V2I transmitter and V2V TX/RX positions.
%
% INPUT
%   d0: highway length
%   d_avg: average inter-vehicle distance, d_avg = 2.5*v
%   Others should be evident from their names
% 
% OUTPUT
%   Flag, 0/1, if the generated vehicles are not enough, Flag = 1
%   vehPos, Nx2 matrix, storing vehicle coordinates, (x, y)
%   indV2I, numV2Ix1, storing indices of V2I transmitters in vehPos
%   indV2Vtx, numV2Vx1, storing indices of V2V transmitters
%   indV2Vrx, numV2Vx1, storing indices of V2V receivers
% 
% Le Liang, Georgia Tech, Aug. 25, 2017


function [ Flag, vehPos, indV2I, indV2Vtx, indV2Vrx ] = genVehLinks(d0, laneWidth, numLane, disBstoHwy, d_avg, numV2I, numV2V)

vehPos = []; % initilizer for vehicle position
indV2I = [];
indV2Vtx = [];
indV2Vrx = [];
Flag = 0;

%% generate all vehicle positions and store in vehPos
for ilane = 1:numLane
    npoints = poissrnd(2*d0/d_avg);
    pproc = (rand(npoints,1)*2-1)*d0; % horizontal coordinates
    pproc2 = [pproc, (disBstoHwy+ilane*laneWidth)*ones(length(pproc), 1)]; % [horizontal vertical] coordinates
    vehPos = [vehPos; pproc2];
end
numVeh = size(vehPos, 1);
if numVeh < numV2I
    Flag = 1; % the generated vehicles are not enough
    return; 
end

%% generate V2I TX and V2V TX/RX positions
indPerm = randperm(numVeh);
indV2I = indPerm(1:numV2I); % randomly pick V2I TXs
indV2I = indV2I(:);
kClosest = numV2V/numV2I;
indV2Vtx = kron(indV2I, ones(kClosest,1));
for ii = 1:numV2I
    vehTmp = vehPos - ones(numVeh,1)*vehPos(indV2I(ii),:); % find the position difference
    [~,ind_ascend] = sort(sum(abs(vehTmp).^2, 2)); % sum along row
    ind_ascend = ind_ascend(:);
    indV2Vrx = [indV2Vrx; ind_ascend(2:kClosest+1)]; % use + 1 to exclude itself when finding the k closest vehicle
end
    
indV2I = indV2I(:);
indV2Vtx = indV2Vtx(:);
indV2Vrx = indV2Vrx(:);


% Use the following commands to show the links graphically. 
% figure;
% plot(vehPos(:,1), vehPos(:,2), 'o');
% hold on
% plot(vehPos(indV2I,1), vehPos(indV2I,2),'k*')
% for ii = 1:numV2V
%     hold on
%     plot(vehPos([indV2Vtx(ii) indV2Vrx(ii)],1), vehPos([indV2Vtx(ii) indV2Vrx(ii)],2),'r-')
% end


    
end

