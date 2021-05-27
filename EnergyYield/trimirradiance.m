% #DESCRIPTION:           This function is used to trim irradiance data
%                         simulated with the Irradiance Module to a
%                         specific wavelength range
%
% #INPUT:                 wavelrange (vector) - wavelength range in nm
%                         I (matrix) - irradiance (size: hours x wavelength)
%                         lambda (vector) - vector of wavelengths which
%                                           represents size(I,2)
%
% #OUTPUT:                I (matrix) - irradiance trimmed to spectral range
%                                      defined by wavelrange
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% -----------------------------------------------------------

function I = trimirradiance(lambda, I, w)
    startindex = find(w==lambda(1));
    stopindex = find(w==lambda(end));
    d = lambda(2)-lambda(1);

    I = I(:,startindex:d:stopindex);
end

