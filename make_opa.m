function [array_ideal, array_nonideal] = make_opa(N, d, amp_sigma, phase_sigma, ...
    angle_steer, lamda, phase_bins)
    %{
    % 1D array
    Inputs:
        N: Number of elements.
        d: Element spacing (m).
        amp_sigma: Amplitude perturbation standard deviation (power, UNITS?).
            Potentially a vector. Normalize to the signal amplitude.
        phase_sigma: Phase perburbation standard deviation (rad).
            Potentially a vector.
        angle_steer: Ideal steering angle (deg);
        lamda: Nominal signal wavelength (m).
        phase_bins: Available phase bins due to driver quantization (rad).
            Values should fall within [0, 2pi]. Note that if you include 0,
            you should include 2pi as well to account for toolbox rounding
            error.
    Outputs:
    %}
    % Constants for future use
    rs = rng(7);
    freq = physconst('lightspeed')/lamda;
    
    % Start with an ideal array
    array_ideal = phased.ULA(N, d);
    
    % Calculate ideal element-to-element phase difference for a target
    % steering angle
    steervec = phased.SteeringVector('SensorArray', array_ideal, ...
                                     'IncludeElementResponse', true);
    sv = steervec(freq, [angle_steer]);
    ind_phase_ideal = angle(sv); % phase setting for each individual element
    ind_phase_ideal = ind_phase_ideal - ind_phase_ideal(1);
    ind_phase_ideal = mod(ind_phase_ideal, 2*pi);
    
    % Find nearest phase match in phase_bins and calculate the phase error
    % for each phased array element (rad)
    ind_phase_error = zeros(size(ind_phase_ideal));
    for i = 1:numel(ind_phase_ideal)
        phase_ideal = ind_phase_ideal(i);
        [~, idx_nearest] = min(abs(phase_ideal-phase_bins));
        phase_nonideal = phase_bins(idx_nearest);
        ind_phase_error(i) = phase_nonideal - phase_ideal;
    end
    
    % Create a phased array which adds amplitude + phase random mismatch and 
    % incorporates the phase error for the same target steering angle as 
    % before
    taper_quant = exp(1i*ind_phase_error).';
    taper_amp = ones(1,N);
    for i=1:numel(amp_sigma)
        taper_amp = taper_amp .* (1 + randn(1,N)*amp_sigma(i));
    end
    
    taper_phase = ones(1,N);
    for i=1:numel(phase_sigma)
        taper_phase = taper_phase .* exp(1i*randn(1,N)*phase_sigma(i));
    end

    array_nonideal = clone(array_ideal);
    array_nonideal.Taper = taper_quant .* taper_amp .* taper_phase;
end