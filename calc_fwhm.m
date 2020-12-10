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
cont_matrix = contour(az_vec, [0,1], pat_vec_norm_concat', [val val]);

% Remove junk from contour
indices = find(cont_matrix(2,:)>0);
cont_matrix(:,indices) = [];
cont_vec = cont_matrix(1,:);

% Find nearest -3dB points relative to peak (making some assumptions here
% about the shape of the output, but they're safe)
[~,k] = mink(abs(cont_vec-az_peak), 2);

if numel(k) < 2
    disp('No FWHM');
    fwhm = Inf;
else
    az_start = cont_vec(min(k));
    az_stop = cont_vec(max(k));
    fwhm = az_stop - az_start;

    if az_start > az_peak
        error('-3dB on left side not included in azimuth sweep');
    end
    if az_stop < az_peak
        error('-3dB point on right side not included in azimuth sweep');
    end
end
end

