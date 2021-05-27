% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function transforms the sun angles phisun
%                         and thetasun into angles described out of the local
%                         rotated coordinate system of the solar cell. All
%                         in- and output angles are in degrees. For more
%                         help see also the functions getillumination() and
%                         EnergyYield().
%
% #INPUT:                 alpha (scalar) - rot angle of solar cell about z
%                         beta (scalar) - rot angle of solar cell about y'
%                         phisun (vector) - the suns' azimuth angles
%                         thetasun (vector) - the suns' solar zenith angles
%
% #OUTPUT:                phisun_rot (vector)
%                         thetasun_rot (vector)
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: eul2quat, quatrotate, round0
%
% -----------------------------------------------------------

function [ phisun_rot, thetasun_rot ] = rotatesunangle( alpha, beta, gamma, phisun, thetasun)

    % first we map the angles to more feasible (mathematical) definitions
    phisun = wrapTo180(phisun);
    thetasun = 90-thetasun;

    % convert to radian & to anglespace used by transformation functions
    phisun = deg2rad(phisun);
    thetasun = deg2rad(thetasun);

    % define euler matrix for rotation: 'ZYX'
    eul = deg2rad([alpha beta gamma]);
    rotm = eul2quat(eul);

    % transform sun coordinates to cartesian
    [sx, sy, sz] = sph2cart(phisun, thetasun, ones(size(phisun)));

    % rotate vectors pointing to sun into rotated solar cell system
    % new coordiants describes sun viewed from rotated local coordinate system
    s = [sx'; sy'; sz']';
    coord_rot = quatrotate(rotm,s);

    % transform back to spherical coordiantes
    [phisun_rot, thetasun_rot, ~] = cart2sph(coord_rot(:,1), coord_rot(:,2), coord_rot(:,3));

    % last transform back to degrees and initial angle definitions
    thetasun_rot = rad2deg(thetasun_rot);
    thetasun_rot = 90-thetasun_rot;
    phisun_rot = rad2deg(phisun_rot);
    phisun_rot = round(wrapTo360(phisun_rot),5);

end

