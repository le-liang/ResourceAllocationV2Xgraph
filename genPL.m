% DESCRIPTION
%   genPL: compute large scale fading (pathloss + shadowing) parameters
%
% INPUT
%   mode = 'V2V' or 'V2I', specify link types of vehcular commun.
%   stdShadow, log-normal shadowing std deviation in dB
%   dist, distance between TX and RX in meters
%   heightTX, transmitter antenna height in meters
%   heightRX, receiver antenna height in meters
%   freq, carrier frequency in GHz
%
% OUTPUT
%   combinedPL, combined large scale fading in dB (pathloss + shadow)
%   
% 
% REFERENCES
%   3rd Generation Partnership Project: Technical Specification Group Radio Access Network: Study LTE-Based V2X Services: (Release 14), 3GPP TR 36.885 V2.0.0, Jun. 2016.
% 
% Le Liang, Georgia Tech, July 10, 2016
% 
% **********************************
% Revised: Le Liang, Georgia Tech, Aug. 25, 2017
%   + add comments and references


function [ combinedPL ] = genPL(linkType, stdShadow, dist, hgtTX, hgtRX, freq)
if strcmp(upper(linkType), 'V2V')
    d_bp = 4*(hgtTX-1)*(hgtRX-1)*freq*10^9/(3*10^8);
    A = 22.7; B = 41.0; C = 20;
    if dist <= 3
        PL = A*log10(3) + B + C*log10(freq/5);
    elseif dist <= d_bp
        PL = A*log10(dist) + B + C*log10(freq/5);
    else
        PL = 40*log10(dist)+9.45-17.3*log10((hgtTX-1)*(hgtRX-1))+2.7*log10(freq/5);
    end
else
    PL = 128.1 + 37.6*log10(sqrt((hgtTX-hgtRX)^2+dist^2)/1000);
end

combinedPL = - (randn(1)*stdShadow + PL);

end