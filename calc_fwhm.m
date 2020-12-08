function [fwhm, az_start, az_stop] = calc_fwhm(az_vec, pat_vec, az_peak)
%%% Inputs:
%%%   az_vec: Azimuth vector (degrees)
%%%   pat_vec: Far-field pattern vector in dB.
%%%   az_peak: Peak about which to calculate the FWHM.
%%% Returns:
%%%   fwhm: Full width half max (in degrees)
%%%   az_start/stop: The angles (in degrees) of intersection for the FWHM

% Normalize the pattern with respect to its peak
pat_peak = interp1(az_vec, pat_vec, az_peak);
pat_vec_norm = pat_vec - pat_peak;

% Find isoline for the -3dB point
val = -3;
pat_vec_norm_concat = [pat_vec_norm, pat_vec_norm];
C = contour(az_vec, [0,1], pat_vec_norm_concat', [val val]);

% Remove junk from contour



end

