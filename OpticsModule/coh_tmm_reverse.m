% ------------- DOKUMENTATION OF THIS FUNCTION -------------
% %DESCRIPTION:           Executes the calculation of the Poynting vector and the
%                         overall scattering matrix for all wavelengths and calculates reflection,
%                         absorption and transmission.
% %INPUT:                 n_array: CRI of the layers (array)
%                         d_list: Thicknesses of the layers (vector)
%                         pol: Polarisation of the field to be computed (string: 'TE' or 'TM')
%                         lam_vac: Wavelength for which the computation shall be conducted (vector)
%                         th_0: Angle of incidendence in degree
% %OUTPUT:                data: detailed overview of the simulation (struct)                       
% %SAVED DATA:            -
% %REQUIRED SUBFUNCTIONS: 
%
% %ADD COMMENTS:          -
% -----------------------------------------------------------

function [R, A, T, th_f] = coh_tmm_reverse(pol, n_array, d_list, th_0, lam_vac)
    % """
    % Reverses the order of the stack then runs coh_tmm.
    % """
    [R, Areverse, T, th_f] = coh_tmm(pol, fliplr(n_array), fliplr(d_list), th_0, lam_vac);
    A = flipud(Areverse);
end