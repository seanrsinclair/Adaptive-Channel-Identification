% ****************************************************************
% RMS phase error from phase noise spectrum
% heavily modified code by A. Apsel
% from code originally written by
% Markus Nentwig, 2011
%
% - generates phase noise spectrum model from single point based on simple approx
% - integrates RMS phase error from a phase noise power spectrum

% ****************************************************************
function [ phi_degRMS] = pn_gen(fcarrier, offset, L1)
    %close all;
    
    % ****************************************************************
    % Phase noise spectrum model
    % --------------------------
    % Notes: 
    % - Generally, a phase noise spectrum tends to follow
    %   |H(f)| = c0 + c1 f^n1 + c2 f^n2 + c3 f^n3 + ...
    %   Therefore, linear interpolation should be performed in dB
    %   on a LOGARITHMIC frequency axis.
    % - Single-/double sideband definition: 
    %   The phase noise model is defined for -inf <= f <= inf
    %   in other words, it contributes the given noise density at
    %   both positive AND negative frequencies. 
    %   Assuming a symmetrical PN spectrum, the model is evaluated 
    %   for |f| => no need to explicitly write out the negative frequency
    %   side.
    % ****************************************************************
    
    % PN model
    % first column: 
    %   frequency offset
    % second column: 
    %   spectral phase noise density, dB relative to carrier in a 1 Hz 
    %   observation bandwidth

    k=1.38e-23;
    T= 290;
    Q = 200;
    fo = 1e5; % 1/f corner
    Psig =1;
    
    %find fitting parameter, F, for Leeson's model from input point
    F=(10^(L1/10))/((2*k*T/Psig)*(1+fo/offset)*(1+(fcarrier./(Q*2*offset))^2));
    
    f_Hz = 10000:100:1e7;
    
    L=10*log10((2*F*k*T/Psig)*(1+fo./f_Hz).*(1+(fcarrier./(Q*2*f_Hz)).^2)); % Leeson's model
    
%     f_Hz= [1 f_Hz];
%     L = [0 L];
    %semilogx(f_Hz, L)

    g_dBc1Hz = L .' -3; %-3dB accounts for 2 sided spectrum of Sphi
    
    % get RMS phase error by integrating the power spectrum
    % (alternative 1)
    phi_degRMS = integratePN_degRMS((f_Hz(2)-f_Hz(1)), g_dBc1Hz);

    % get RMS phase error based on a test signal
    % (alternative 2)
%     n = 2 ^ 20;
%     deltaF_Hz = 2;    
%     e2_degRMS = simulatePN_degRMS(f_Hz, g_dBc1Hz, n, deltaF_Hz)
end

% ****************************************************************
% Integrate the RMS phase error from the power spectrum
% ****************************************************************
function r = integratePN_degRMS(f_Hz, g1_dBc1Hz)
    1;
    % integration step size
    deltaF_Hz = f_Hz(1);
    
    % list of integration frequencies
%     f_Hz = deltaF_Hz:deltaF_Hz:5e6;
    
    % interpolate spectrum on logarithmic frequency, dB scale
    % unit is dBc in 1 Hz, relative to a unity carrier
%     fr_dB = interp1(log(f_Hz(1)+eps), g1_dBc1Hz, log(f_Hz+eps), 'linear');
    
    % convert to power in 1 Hz, relative to a unity carrier
    fr_pwr = 10 .^ (g1_dBc1Hz/10);

    % scale to integration step size
    % unit: power in deltaF_Hz, relative to unity carrier
    fr_pwr = fr_pwr * deltaF_Hz;

    % evaluation frequencies are positive only
    % phase noise is two-sided 
    % (note the IEEE definition: one-half the double sideband PSD)    
    fr_pwr = fr_pwr * 2;

    % sum up relative power over all frequencies
    pow_relToCarrier = sum(fr_pwr);

    % convert the phase noise power to an RMS magnitude, relative to the carrier
    pnMagnitude = sqrt(pow_relToCarrier);

    % convert from radians to degrees
    r = pnMagnitude * 180 / pi;
end


