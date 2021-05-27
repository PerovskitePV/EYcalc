% ------------- DOCUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function calculates the rotation angle alpha
%                         and tilt angle beta from a constant tilt angle and
%                         the azimuth and elevation of sun.
%
%                         ATTENTION: status of this function is "testing"
%
%
% #INPUT:                 tilt (scalar) - const tilt angle, cell facing south
%                         phi (scalar) - azimuth of sun
%                         theta (scalar) - elevation of sun
%
% #OUTPUT:                alpha (scalar) - rotation angle of cell
%                         beta (scalar) - tilt angle of cell
%                         gamma (scalar) - 3rd euler angle of cell rotation
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% -----------------------------------------------------------

function [ alpha, beta, gamma ] = seesawtilt( tilt, rotation, phisun0, thetasun0 )

    alpha = rotation;
    beta = tilt;
    
    % Eueler matrix
    eul = deg2rad([alpha beta 0]);
    rotm = eul2quat(eul,'ZYX');

    % rotate normal of cell with quadrotate. Be carefull: quaternions use
    % right handed coordinate systems!
    v = quatrotate(rotm,[0 0 1]);
    u = quatrotate(rotm,[0 1 0]);

    % Definition of sun position
    x = [cosd(phisun0)*sind(thetasun0), -sind(phisun0)*sind(thetasun0), cosd(thetasun0)];

    % Orthogonal projection into rotation plane
    P = x*u'/(u*u')*u + x*v'/(v*v')*v;
    
    % seesaw rotation angle
    gamma = acosd(P*v'/(norm(P)*norm(v)));
    if phisun0 >= 180
        gamma = -gamma;
    end

end
