% ------------- DOCUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           Calculates direction of transmitted ray through
%                         the facet of a pyramidal structure for all wavelengths.
% #INPUT:                 v: Vector of incident light ray (vector)
%                         n: Surface normal of facet (vector)
%                         n1,n2: Refractive indices of starting and end layer (array)
% #OUTPUT:                vtrans: Vector of transmitted ray (vector)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          -
% -----------------------------------------------------------

function vtrans = snell3D(v,n,n1,n2)
    %% refract ray with velocity v at plane with normal vector n separating regions with refractive indices n1 and n2
    % normalize input vectors and ensure n is pointing upwards (against
    % direction of v)
    v = v/norm(v);
    n = -n/norm(n)*sign(dot(n,v));
    
    % calculate direction of refracted ray
    vtrans = n1./n2.*repmat((cross(n,cross(-n,v)))',1,length(n1))-n'*sqrt(1-(n1./n2).^2.*repmat(norm(cross(n,v))^2,1,length(n1)));
       
    % normalize
    vtrans = vtrans./vecnorm(vtrans);
end