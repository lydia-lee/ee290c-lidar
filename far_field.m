clear all;
close all;
clc;

%% Create array
% --------- BEGIN EDIT ---------
N = 1000;  % number of elements
d = 7.4e-6; %[1.55e-6/2, 1.55e-6/2]; % element spacing (m)
sigma_taper = (d*N/2)/3; % standard deviation of taper shape (m)
mu_taper = 0; % location of center of the normal beam (m)
% ---------- END EDIT ----------

array = phased.ULA(N, d);
% viewArray(array,'Title','Uniform Linear Array (ULA)')

% Generating taper
% array_tapered = clone(array);
% array_tapered.Taper = gausswin(N);
% array = array_tapered;

%% Simulate antenna behavior
% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % wavelength (m)
theta_steer = [6]; % steering angle azimuth;elevation (deg)

az = 5:0.00001:7; % azimuth angles (degrees)
el = 0; %-90:1:90; % elevation angles (degrees)
coord = 'rectangular';
plot_quant = 'powerdb';
% ---------- ENDEDIT ----------

psi_shift = (sin(theta_steer*pi/180)*2*pi.*d / lamda) * 180/pi; % phase shift in each adjacent element (degrees)
Lamda = psi_shift*1/360.*d; % reset distance (m)
freq = physconst('lightspeed')/lamda;
steervec = phased.SteeringVector('SensorArray', array, ...
                                 'PropagationSpeed', physconst('lightspeed'), ...
                                 'IncludeElementResponse', false, ...
                                 'NumPhaseShifterBits', 4);

sv = steervec(freq, theta_steer);
[pat_vec, az_vec, el_vec] = pattern(array, freq, az, el, 'PropagationSpeed', physconst('lightspeed'), ...
                             'CoordinateSystem', coord, ...
                             'Type', plot_quant, ...
                             'Weights', sv);
if isscalar(el_vec)
    
    for i=1:size(pat_vec, 2)
        subplot(size(pat_vec, 2), 1, i);
        pat = pat_vec(1:end, i);
        plot(az_vec, pat, 'DisplayName', sprintf('%0.1f deg', theta_steer(i)));
        
        [pks_unsort, idx_pks_unsort] = findpeaks(pat);
        
        % Getting the main lobe (assuming steering angle is correct)
        [scratch, idx_real_main] = min(abs(az_vec-theta_steer(i)));
        az_main = az_vec(idx_real_main);
        
        % Finding next largest sidelobe
        [pks, idx_idx] = sort(pks_unsort, 'descend');
        idx_main = idx_idx(1);
        idx_side = idx_idx(2);
        az_side = az_vec(idx_pks_unsort(idx_side));
        % text(az_main, pks(1), sprintf('%0.1f dB @ %0.1f deg', pks(1), az_main));
%         text(az_side, pks(2), sprintf('%0.1f dB @ %0.1f deg', pks(2), az_side));
        
        % Finding nearest sidelobe
        pks = pks(2:end);
        idx_idx = idx_idx(2:end);
        idx_idx_shifted = idx_idx-idx_main;
        idx_nearest = min(abs(idx_idx_shifted)) + idx_main;
        az_nearest = az_vec(idx_pks_unsort(idx_nearest));
        pk_nearest = pat(idx_pks_unsort(idx_nearest));
%         text(az_nearest, pk_nearest, sprintf('%0.1f dB @ %0.1f^o', pk_nearest, az_nearest))
        xlabel("Azimuth (^o)");
        ylabel("Power (dB)");
        title(sprintf("%0.1f^o Steering, %0.0f Elements", theta_steer(i), N));
        ylim([-60, 0.1]);
        
        % Finding the FWHM (1/e) of the intended main lobe
        val = -3; % 10*log10(1/exp(1));
        idx_pat = find(diff(sign(pat-val)))+1;
        fwhm = abs(az_vec(idx_pat(2)) - az_vec(idx_pat(1)));
        fwhm_rad = fwhm*pi/180;
        text(az_main, val, sprintf("FWHM %0.2e rad", fwhm_rad));
        hold on;
        plot(az_vec, val*ones(size(az_vec)));
end

else
    mesh(az_vec, el_vec, pat_vec);
    xlabel("Azimuth (deg)");
    ylabel("Elevation (deg)");
    zlabel("Power (dB)");
end