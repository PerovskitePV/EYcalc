% ------------- DOCUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           Calculates absorption based on refractive index and thicknesses of a
%                         layer for all angles of incidence. t needs to be a scalar, n needs to be
%                         a vector.
% #INPUT:                 n: Refractive indices of incoherent layer (array)
%                         t: Layer thicknesses of effective interface (array)
%                         lambda: Wavelengths (vector)
% #OUTPUT:                abs: Absorption of incoherent layer for all wavelengths and all angles of incidence(array)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          
%                         
% -----------------------------------------------------------

function abs = Lambert(n,t,lambda)
    abs = (1-exp(-(1./cosd(0:89))'*(4*pi./lambda).*imag(n)*t));
end
