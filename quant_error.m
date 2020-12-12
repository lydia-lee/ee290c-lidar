clear all;
close all;
clc;

% ---- BEGIN HERE ----
N_sig = 2; % number of bits used to quantize the controlling signal (e.g. voltage, current, etc.)
N_phase = 2; % number of bits ideally used to quantize the phase shift (e.g. 1 bit -> 0 and pi, 2 bits = 0, pi/2, pi, 3pi/2)
phase_sig_pow = 2; % exponent of signal-to-phase relationship, e.g. 2 --> phase is proportional to (sig)^2
sig_2pi = 1; % signal value which corresponds to 2pi phase shift
FSR = 1; % FSR of the DAC
sig_type = 'Signal'; % Used for a plotting label
% ----- END HERE -----

LSB_phase = 2*pi/2^N_phase;
LSB_phase_deg = LSB_phase * 180/pi;
LSB_sig = FSR/2^N_sig;

scale = sig_2pi/((2*pi)^(1/phase_sig_pow)); % to match the 2pi point

phase_quant_ideal = [0:2^N_phase-1] * LSB_phase; % ideal quantized phase
sig_ideal = phase_quant_ideal.^(1/phase_sig_pow) * scale; % ideal signal value corresponding to each phase value
sig_quant = [0:2^N_sig-1] * LSB_sig; % available quantized signal values
phase_available = (sig_quant./scale).^(phase_sig_pow);

% Finding the closest approximation possible for each ideal signal value
sig_real = zeros(size(sig_ideal));
for i=1:length(sig_ideal)
    [scratch, idx_closest] = min(abs(sig_ideal(i)-sig_quant));
    sig_closest = sig_quant(idx_closest);
    sig_real(i) = sig_closest;
end

% Back-calculating the phase shift error this introduces
phase_quant_real = (sig_real./scale).^(phase_sig_pow);

%%% Plotting results
% Append final value to the end for a cleaner plot
sig_ideal = cat(2, sig_ideal, sig_ideal(end));
sig_real = cat(2, sig_real, sig_real(end));
phase_quant_ideal = cat(2, phase_quant_ideal, phase_quant_ideal(end));
phase_quant_real = cat(2, phase_quant_real, phase_quant_real(end));

figure;

subplot(3, 1, 1);
    stairs(sig_ideal, 'Linewidth', 2);
    hold on;
    stairs(sig_real, 'Linewidth', 2);
    legend('Ideal', 'Real');
    ylabel(sig_type);
    xlim([1,length(sig_ideal)]);
    
[t, s] = title({"\phi_{step} = " + sprintf("%0.0d^o", LSB_phase_deg), ...
         sprintf("N_{bits} = %0.0d", N_sig) + ", " + ...
         sprintf("V_{FSR} = %0.1fV_{reset}", FSR/sig_2pi)});
t.FontSize = 14;    

subplot(3, 1, 2);
    stairs(phase_quant_ideal/pi, 'Linewidth', 2);
    hold on;
    stairs(phase_quant_real/pi, 'Linewidth', 2);
    legend('Ideal', 'Real');
    ylabel('Phase (\pi rad)');
    xlim([1,length(phase_quant_ideal)]);
    
subplot(3, 1, 3);
    stairs((phase_quant_ideal - phase_quant_real)/pi, 'Linewidth', 2);
    ylabel('Phase Error (\pi rad)');
    xlim([1,length(phase_quant_ideal)]);