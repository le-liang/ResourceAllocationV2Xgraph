% DESCRIPTION
%   compute the optimal V2I and V2V power allocations given one V2I link
%   and multiple V2V links sharing spectrum
%
% INPUT
%   alpha_m, numK x 1, interfering channel from V2I to numK V2Vs
%   alpha_kk, numK x numK, mutual V2V interference channels
%   alpha_mB, scalar, V2I channel
%   alpha_kB, numK x 1, interfering channel from numK V2Vs to V2I link
%   gamma0, scalar, V2V slow SINR threshold
%   sig2, noise power
%   Pd_max, V2V max power
%   Pc_max, V2I max power
%
% OUTPUT
%   rate, maximum rate achieved with the optimal power allocation, -1 if
%   the problem is infeasible. 
%   Pc_opt, optimal power of V2I
%   Pd_opt, numK x 1, optimal power of V2Vs
%
% Le Liang, Georgia Tech, Aug. 30, 2017



function [rate, Pc_opt, Pd_opt] = calOptPower(alpha_mk, alpha_kk, alpha_mB, alpha_kB, gamma0, sig2, Pd_max, Pc_max)

numK = length(alpha_mk); % # of mutually interfering V2Vs

phi = alpha_kk';
phi = -gamma0*phi;
for ii = 1:numK
    phi(ii,ii) = alpha_kk(ii,ii);
end

phi_inv = inv(phi);
Pc_cand = (Pd_max-gamma0*sig2*sum(phi_inv,2))./(gamma0*phi_inv*alpha_mk);

Pc_opt = min([Pc_max; Pc_cand]);
if Pc_opt <= 0 % infeasible
    rate = -1;
    Pd_opt = 0;
    return;
end
Pd_opt = phi_inv*(gamma0*(Pc_opt*alpha_mk + sig2));
if sum(Pd_opt<=0) >= 1 % infeasible
    rate = -1;
    return;
end

signal = Pc_opt*alpha_mB;
interf = Pd_opt'*alpha_kB;% interference
rate = log2(1 + signal/(sig2 + interf));
end






