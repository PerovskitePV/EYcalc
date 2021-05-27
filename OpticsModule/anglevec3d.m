% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function calculates the anlge between two vectors in 3d
%
% #INPUT:                 a - vector 1 as [x y z] or matrix
%                         b - vector 2 as [x y z] or matrix
%                         c - boolean, if true take absolute value of dot product
%
% #OUTPUT:                angle as scalar or vector
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% -----------------------------------------------------------

function angle = anglevec3d(a,b,f)  
    % format input
    if size(a,1) ~= 3
       a = shiftdim(a,1); 
    end
    if size(b,1) ~= 3
       b = shiftdim(b,1);
    end
    
    % bring to same size / dimensions
    a = repmat(a,1,size(b,2));
    b = repmat(b,1,size(a,2));
    
    if f==true
        angle = floor(atan2d( vecnorm(cross(a,b)), abs(dot(a,b)) ));  
    else
        angle = floor(atan2d( vecnorm(cross(a,b)), dot(a,b) )); 
    end
    
    % Other possible formulations, leading to worse/false results:    
    % alphaTransIdx = floor(round(abs(atand(vtrans(1,:)./vtrans(3,:))),10))+1; % angle of transmission vector
    % alphaTransIdx = floor(acosd( dot(vtrans,repmat([0 0 -1]',1,size(vtrans,2))) ))+1;
    
end

