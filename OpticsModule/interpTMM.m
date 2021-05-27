% ------------- DOCUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           Calculates absorption based on refractive index and thicknesses of a
%                         layer for all angles of incidence. t needs to be a scalar, n needs to be
%                         a vector.
% #INPUT:                 Ref: Reflection calculated by TMM (array)
%                         Abs: Absorption calculated by TMM (tensor)
%                         Trans: Transmission calculated by TMM (array)
%                         n1,n2: Refractive indices of starting and end layer (vector)
%                         thickness: Thicknesses of TMM layers (vector)
%                         TMMAngleRange: AngleRange of TMM calculation (vector)
%                         lambda: Wavelengths (vector)
% #OUTPUT:                abs: Absorption of incoherent layer for all wavelengths and all angles of incidence(array)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: transAngles
%
% #ADD COMMENTS:          
%                         
% -----------------------------------------------------------

function [RefInterp,AbsInterp,TransInterp,TransAngles] = interpTMM(Ref,Abs,Trans,n1,n2,thickness,TMMAngleRange,lambda)
% Function interpTMM
% This function interpolates the TMM results to 1 deg of incidence and also
% considers total internal reflection
    AngleRange = 0:89;
    TransAngles= floor(transAngles(n1,n2));
    RefInterp = NaN(length(AngleRange),length(lambda));
    AbsInterp = NaN(length(thickness)-2,length(AngleRange),length(lambda));
    TransInterp = NaN(length(AngleRange),length(lambda));
    RefInterp(TMMAngleRange+1,:) = Ref;
    AbsInterp(:,TMMAngleRange+1,:) = Abs;
    TransInterp(TMMAngleRange+1,:) = Trans;
    RefInterp = fillmissing(RefInterp,'linear');
    for i=1:length(thickness)-2
        AbsInterp(i,:,:) = fillmissing(reshape(AbsInterp(i,:,:),length(AngleRange),length(lambda)),'linear');
    end
    TransInterp = fillmissing(TransInterp,'linear');
end