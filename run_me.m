clear all;
close all;
clc;
% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % meters
N = 2000;
d = 4.5e-6; % meters
amp_sigma = [0, 0]; % normalize to 1
phase_sigma = [0, 0]; % rad
angle_steer = 10; % deg
phase_bins = [0,0.0245,0.0982,0.2209,0.3927,0.6136,0.8836,1.2026,1.5708,1.9880,2.4544,2.9698,3.5343,4.1479,4.8106,5.5223,2*pi];
% [0,0.3927,1.5708,3.5343]; % [0,0.0245,0.0982,0.2209,0.3927,0.6136,0.8836,1.2026,1.5708,1.9880,2.4544,2.9698,3.5343,4.1479,4.8106,5.5223,2*pi]; % rad
angle_plot =  0:0.001:angle_steer*1.1; % azimuth angles (degrees)
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

d_ideal = directivity(array_ideal, freq, angle_steer, ...
                        'Weights', sv_ideal);
d_nonideal = directivity(array_nonideal, freq, angle_steer, ...
                         'Weights', sv_ideal);

disp("Directivity Ideal/Nonideal: " + d_ideal + "/" + d_nonideal + " dB");                     
                     
fwhm_ideal = calc_fwhm(az_ideal_vec, pat_ideal_vec, angle_steer);
fwhm_nonideal = calc_fwhm(az_nonideal_vec, pat_nonideal_vec, angle_steer);
disp("FWHM Ideal/Nonideal: " + fwhm_ideal*pi/180 + "/" + fwhm_nonideal*pi/180 + " radians");

% Plot array outputs over one another
% plot(az_ideal_vec, pat_ideal_vec, 'DisplayName', sprintf('%0.1f deg', angle_steer));
% figure;
plot(az_nonideal_vec, pat_nonideal_vec, 'DisplayName', sprintf('%0.1f deg', angle_steer));
xlim([min(angle_plot), max(angle_plot)])
% legend('Ideal', 'Nonideal');

%%
% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % meters
N = 20000;
d_vec = [1.0:0.1:4.5] * 1e-6; % meters - 9.0
amp_sigma = [0, 0]; % normalize to 1
phase_sigma = [0, 0]; % rad
angle_steer = 10; % deg
phase_bins = [0, pi, 2*pi]; % rad
fwhm_max = 200e-6; % rad
angle_plot =  9.8:0.00001:10.2; % azimuth angles (degrees)
% ---------- END EDIT ----------
c = physconst('lightspeed');
freq = c/lamda;
N_vec = zeros(size(d_vec));

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
        
        % Break if you get too high
        if N_mid > 30000
            disp("break @ d=" + d);
            break
        end
        
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
                N_low = N_mid;
                N_high = N_high + N_mid;
            else
                range_set = true;
            end

        % Binary search after we've overshot the spec
        else
            if N_high <= N_mid || N_low >= N_mid
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

disp(d_vec);
disp(N_vec);

%%
% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % meters
N = 20000;
d_vec = [1.0:0.1:4.5] * 1e-6; % meters - 9.0
amp_sigma = [0, 0]; % normalize to 1
phase_sigma = [0, 0]; % rad
angle_steer = 10; % deg
phase_bins = [0,0.0245,0.0982,0.2209,0.3927,0.6136,0.8836,1.2026,1.5708,1.9880,2.4544,2.9698,3.5343,4.1479,4.8106,5.5223,2*pi]; % rad
fwhm_max = 100e-6; % rad
angle_plot =  9.8:0.00001:10.2; % azimuth angles (degrees)
% ---------- END EDIT ----------
c = physconst('lightspeed');
freq = c/lamda;
N_vec = zeros(size(d_vec));

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
        
        % Break if you get too high
        if N_mid > 30000
            disp("break @ d=" + d);
            break
        end
        
        [array_ideal, array_nonideal] = make_opa(N_mid, d, amp_sigma, phase_sigma, ...
            angle_steer, lamda, phase_bins);

        % Crate steering vector for ideal array
        steervec_ideal = phased.SteeringVector('SensorArray', array_ideal, ...
            'PropagationSpeed', c, ...
            'IncludeElementResponse', true, ...
            'NumPhaseShifterBits', 0);
        sv_ideal = steervec_ideal(freq, angle_steer);

        % Apply steering vector to ideal array
        [pat_nonideal_vec, az_nonideal_vec, ~] = pattern(array_nonideal, freq, angle_plot, 0, ...
            'PropagationSpeed', c, ...
            'CoordinateSystem', 'rectangular', ...
            'Type', 'powerdb', ...
            'Weights', sv_ideal);

        fwhm = calc_fwhm(az_nonideal_vec, pat_nonideal_vec, angle_steer) * pi/180;

        % Increase N_high and N_low until we've overshot the spec
        if ~range_set
            if fwhm > fwhm_max
                N_low = N_mid;
                N_high = N_high + N_mid;
            else
                range_set = true;
            end

        % Binary search after we've overshot the spec
        else
            if N_high <= N_mid || N_low >= N_mid
                disp("d/N: " + d + "/" + N_mid);
                N_vec(i) = N_mid;
                done = true;
                figure;
                plot(az_nonideal_vec, pat_nonideal_vec);
            elseif fwhm < fwhm_max
                N_high = N_mid;
            else
                N_low = N_mid;
            end

        end
    end
end

disp(d_vec);
disp(N_vec);

%%
clear all;
close all;
clc;
% --------- BEGIN EDIT ---------
lamda = 1.55e-6; % meters
N = 1850; % number of array elements
d = 7.4e-6; % meters
amp_sigma = [0, 0]; % normalize to 1
phase_sigma = [0, 0]; % rad
angle_steer = 6; % deg

N_sig_vec = 1:6;
phase_sig_pow_vec = [1,2,4];
overshoot_vec = 0:0.1:0.5;
% ---------- END EDIT ----------
c = physconst('lightspeed');
freq = c/lamda;
angle_plot =  0:0.001:angle_steer*1.1; % azimuth angles (degrees)

for N_sig = N_sig_vec
    for phase_sig_pow = phase_sig_pow_vec
        for overshoot = overshoot_vec
            phase_bins = quant_phase_values(N_sig, phase_sig_pow, 1, overshoot);
            disp("N_sig/phase_sig_pow/overshoot: " + N_sig + "/" + phase_sig_pow + "/" + overshoot);

            % Create ideal array and its counterpart with noise (includes quantized
            % control, etc.)
            [array_ideal, array_nonideal] = make_opa(N, d, amp_sigma, phase_sigma, ...
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

            % Apply same steering vector to noisy array
            [pat_nonideal_vec, az_nonideal_vec, ~] = pattern(array_nonideal, freq, angle_plot, 0, ...
                'PropagationSpeed', c, ...
                'CoordinateSystem', 'rectangular', ...
                'Type', 'powerdb', ...
                'Weights', sv_ideal);
    
%             pks_nonideal = findpeaks(pat_nonideal_vec);
%             sidelobe_nonideal = min(maxk(pks_nonideal, 2));
%             disp("Normalized Sidelobe: " + sidelobe_nonideal + " dB");
            
            d_ideal = directivity(array_ideal, freq, angle_steer, ...
                                    'Weights', sv_ideal);
            d_nonideal = directivity(array_nonideal, freq, angle_steer, ...
                                     'Weights', sv_ideal);

            disp("Directivity Ideal/Nonideal: " + d_ideal + "/" + d_nonideal + " dB");

%             fwhm_ideal = calc_fwhm(az_ideal_vec, pat_ideal_vec, angle_steer);
%             fwhm_nonideal = calc_fwhm(az_nonideal_vec, pat_nonideal_vec, angle_steer);
%             disp("FWHM Ideal/Nonideal: " + fwhm_ideal*pi/180 + "/" + fwhm_nonideal*pi/180 + " radians");

            % Plot array outputs over one another
            % plot(az_ideal_vec, pat_ideal_vec, 'DisplayName', sprintf('%0.1f deg', angle_steer));
%             figure;
%             plot(az_nonideal_vec, pat_nonideal_vec, 'DisplayName', sprintf('%0.1f deg', angle_steer));
%             xlim([min(angle_plot), max(angle_plot)]);
%             title({"N_{sig} = " + N_sig, ...
%                     "phase\_sig\_pow = " + phase_sig_pow, ...
%                     "overshoot = " + overshoot});
            % legend('Ideal', 'Nonideal');
        end
    end
end