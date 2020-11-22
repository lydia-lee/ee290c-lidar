clear all;
close all;
clc;

% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % meters
N = 1000;
d = 7.4e-6; % meters
amp_sigma = [0]; % normalize to 1
phase_sigma = [0]; % rad
angle_steer = 6; % deg
phase_bins = [0, pi/2]; % rad
angle_plot =  5:0.00001:7; % azimuth angles (degrees)
% ---------- END EDIT ----------
c = physconst('lightspeed');
freq = c/lamda;

[array_ideal, array_nonideal] = make_opa(N, d, amp_sigma, phase_sigma, ...
    angle_steer, lamda, phase_bins);

steervec = phased.SteeringVector('SensorArray', array_nonideal, ...
    'PropagationSpeed', c, ...
    'IncludeElementResponse', true, ...
    'NumPhaseShifterBits', 0);

sv = steervec(freq, angle_steer);
[pat_vec, az_vec, el_vec] = pattern(array_nonideal, freq, angle_plot, 0, ...
    'PropagationSpeed', c, ...
    'CoordinateSystem', 'rectangular', ...
    'Type', 'powerdb', ...
    'Weights', sv);

for i=1:size(pat_vec, 2)
    subplot(size(pat_vec, 2), 1, i);
    pat = pat_vec(1:end, i);
    plot(az_vec, pat, 'DisplayName', sprintf('%0.1f deg', angle_steer(i)));
end