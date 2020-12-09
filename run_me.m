clear all;
close all;
clc;
%%
% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % meters
N = 3000;
d = 4.5e-6; % meters
amp_sigma = [0, 0]; % normalize to 1
phase_sigma = [0, 0]; % rad
angle_steer = 10; % deg
phase_bins = [0, pi, 2*pi]; % rad
angle_plot =  9:0.0001:11; % azimuth angles (degrees)
% ---------- END EDIT ----------
c = physconst('lightspeed');
freq = c/lamda;

% Create ideal array and its counterpart with noise (includes quantized
% control, etc.)
[array_ideal, array_nonideal] = make_opa(N, d, amp_sigma, phase_sigma, ...
    angle_steer, lamda, phase_bins);

% Crate steering vector for ideal array
steervec_ideal = phased.SteeringVector('SensorArray', array_ideal, ...
    'PropagationSpeed', c, ...
    'IncludeElementResponse', true, ...
    'NumPhaseShifterBits', 3);
sv_ideal = steervec_ideal(freq, angle_steer);

% Apply steering vector to ideal array
[pat_ideal_vec, az_ideal_vec, ~] = pattern(array_ideal, freq, angle_plot, 0, ...
    'PropagationSpeed', c, ...
    'CoordinateSystem', 'rectangular', ...
    'Type', 'powerdb', ...
    'Weights', sv_ideal);

% Apply same steering vector to noisy array
[pat_nonideal_vec, az_nonideal_vec, ~] = pattern(array_nonideal, freq, angle_plot, 0, ...
    'PropagationSpeed', c, ...
    'CoordinateSystem', 'rectangular', ...
    'Type', 'powerdb', ...
    'Weights', sv_ideal);

fwhm_ideal = calc_fwhm(az_ideal_vec, pat_ideal_vec, angle_steer);
disp("FWHM: " + fwhm_ideal*pi/180 + " radians");

% Plot array outputs over one another
plot(az_ideal_vec, pat_ideal_vec, 'DisplayName', sprintf('%0.1f deg', angle_steer));
% hold on;
% plot(az_nonideal_vec, pat_nonideal_vec, 'DisplayName', sprintf('%0.1f deg', angle_steer));
% legend('Ideal', 'Nonideal');

%%
% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % meters
N = 500;
d_vec = [1.0:0.5:20] * 1e-6; % meters
amp_sigma = [0, 0]; % normalize to 1
phase_sigma = [0, 0]; % rad
angle_steer = 10; % deg
phase_bins = [0, pi, 2*pi]; % rad
fwhm_max = 100e-6; % rad
angle_plot =  0:0.0001:11; % azimuth angles (degrees)
% ---------- END EDIT ----------
c = physconst('lightspeed');
freq = c/lamda;
N_vec = ones(size(d_vec));

for i = 1:numel(d_vec)
    d = d_vec(i);
    N_low = 0;
    N_high = N;
    range_set = false;
    done = false;
    % Create ideal array and its counterpart with noise (includes quantized
    % control, etc.)
    while ~done
        N_mid = round((N_high+N_low)/2);
        [array_ideal, ~] = make_opa(N_mid, d, amp_sigma, phase_sigma, ...
            angle_steer, lamda, phase_bins);

        % Crate steering vector for ideal array
        steervec_ideal = phased.SteeringVector('SensorArray', array_ideal, ...
            'PropagationSpeed', c, ...
            'IncludeElementResponse', true, ...
            'NumPhaseShifterBits', 0);
        sv_ideal = steervec_ideal(freq, angle_steer);

        % Apply steering vector to ideal array
        [pat_ideal_vec, az_ideal_vec, ~] = pattern(array_ideal, freq, angle_plot, 0, ...
            'PropagationSpeed', c, ...
            'CoordinateSystem', 'rectangular', ...
            'Type', 'powerdb', ...
            'Weights', sv_ideal);

        fwhm = calc_fwhm(az_ideal_vec, pat_ideal_vec, angle_steer) * pi/180;

        % Increase N_high and N_low until we've overshot the spec
        if ~range_set
            if fwhm > fwhm_max
                N_low = N_high;
                N_high = N_high * 2;
            else
                range_set = true;
            end

        % Binary search after we've overshot the spec
        else
            if fwhm == fwhm_max || N_high <= N_mid || N_low >= N_mid
                disp("d/N: " + d + "/" + N_mid);
                N_vec(i) = N_mid;
                done = true;
            elseif fwhm < fwhm_max
                N_high = N_mid;
            else
                N_low = N_mid;
            end

        end
    end
end