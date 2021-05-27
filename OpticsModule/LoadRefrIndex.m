% ------------- DOKUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           This function returns the complex index of refraction 
%                         spectra, ntotal, for the material called 'name' for 
%                         each wavelength value in the wavelength vector 'wavelengths'.  
%                         The material must be present in the index of refraction
%                         library 'Index_of_Refraction_library.xls'.  The program uses linear
%                         interpolation/extrapolation to determine the index of refraction for
%                         wavelengths not listed in the library.
% #INPUT:                 AngleRange: Range of incident angle (vector)
%                         Morphology: Surface morphology (string)
% #OUTPUT:                ntotal: Complex refractive index for all
%                         materials and all wavelengths (tensor)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          -
% -----------------------------------------------------------
function ntotal = LoadRefrIndex(name,wavelengths,IndRefr,IndRefr_names)
    try 
        % Load index of refraction data in spread sheet, will crash if misspelled
        file_wavelengths=IndRefr(:,strcmp('Wavelength (nm)',IndRefr_names));
        n=IndRefr(:,strcmp(strcat(name,'_n'),IndRefr_names));
        k=IndRefr(:,strcmp(strcat(name,'_k'),IndRefr_names));

        % Interpolate/Extrapolate data linearly to desired wavelengths
        n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
        k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

        %Return interpolated complex index of refraction data
        ntotal = n_interp+1i*k_interp;
    catch
        error(['Material: ',name,' not found in RefractiveIndexLibrary']);
    end
end