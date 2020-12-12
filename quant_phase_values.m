function [phase_available] = quant_phase_values(N_sig, phase_sig_pow, sig_2pi, overshoot)
%%% Inputs:
%%%     N_sig: Number of bits for the controlling DAC.
%%%     phase_sig_pow: Exponent of signal-to-phase relationship, e.g. 1 -->
%%%         linear, 2 --> phase proportional to sig^2, etc.
%%%     sig_2pi: Signal quantity which corresponds to a 2pi phase shift.
%%%     overshoot: Fraction by which the controller max output overshoots
%%%         sig_2pi. e.g. overshoot=0.1 --> DAC FSR is 
%%%         1.1x what's required to reach 2pi modulation.
%%% Returns:
%%%     phase_available: The phase values available given the input
%%%     parameters.
%%% Example:
%%%     N_sig = 2, phase_sig_pow = 1, sig_2pi = 1, overshoot = 0.
%%%     returns [0, pi/2, pi, 3pi/2]
FSR = sig_2pi * (1+overshoot);
LSB_sig = FSR/2^N_sig;

scale = sig_2pi/((2*pi)^(1/phase_sig_pow)); % to match the 2pi point

sig_quant = [0:2^N_sig-1] * LSB_sig; % available quantized signal values
phase_available = (sig_quant./scale).^(phase_sig_pow);

end