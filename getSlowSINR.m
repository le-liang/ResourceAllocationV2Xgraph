% DESCRIPTION
%   Get equivalent SINR expression in terms of slow CSI according to 
%   Eq. (3) in [R1] while satisfying outage probability
%
% INPUT
%   gamma_tar, instantaneous SINR threshold target, in dB
%   outage, acceptable outage probability
%
% OUTPT
%   gamma, SINR expressed in slow CSI, in dB
%
% REFERENCE
%   [R1] W. Sun, et al "Cluster-based radio resource management for
%   D2D-supported safety-critical V2X communications," IEEE Trans. Wireless
%   Commun., vol. 15, no. 4, pp. 2756--2769, Apr. 2016.
%
% Le Liang, Georgia Tech, Aug. 25, 2017
%
%
% Note: 
%   According to [R1], gamma_tar = gamma/(-log(1-outage)). Thus the
%   function can be simplified. 
%

function [gamma] = getSlowSINR(gamma_tar, outage)

gamma_tar = 10^(gamma_tar/10);
gamma = gamma_tar/(-log(1-outage));
gamma = 10*log10(gamma);

% num = 1e7;
% hi = (randn(num,1)+1j*randn(num,1))/sqrt(2);
% 
% x = zeros(num,1);
% for cnt = 1:num
%     x(cnt) = abs(hi(cnt))^2;
% end
% 
% [y] = sort(x); % sort x in ascending order
% 
% ind = floor(outage*num);
% 
% gamma = gamma_tar/y(ind);
% gamma = 10*log10(gamma);



