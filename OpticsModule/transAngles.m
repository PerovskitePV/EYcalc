% ------------- DOCUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           This function calculates the transmission angle at an interface 
%                         between materials with complex dielectric constant n1 and n2, for all 
%                         incident angles theta_i
% #INPUT:                 n1,n2: Refractive indices of starting and end layer (vector)
% #OUTPUT:                theta: Transmission angles (vector)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          
%                         
% -----------------------------------------------------------

function theta = transAngles(n1,n2)
    theta_i = 0:89;
    theta = asind(repmat(real(n1)./real(n2),[length(theta_i) 1]).*repmat(sind(theta_i).',[1 length(n1)]));
end