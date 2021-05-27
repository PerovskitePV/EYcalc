% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function calculates a matrix GI with ones
%                         and zeros. Ones are at corresponding theta and
%                         phi (local spherical coordinates) pairs where
%                         sun-rays can hit the solar cell. Zeros define
%                         angle combinations from where no light can hit
%                         the solar cell. This can be due to some rotation.
%                         The solar cell rotation angles are defined by
%                         alpha (about azimuth) and then about beta
%                         (rotation about new local y axis). See
%                         EnergyYield for more help.
%
% #INPUT:                 alpha (scalar) - rot angle of solar cell about z
%                         beta (scalar) - rot angle of solar cell about y'
%
% #OUTPUT:                GI (matrix)
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: eul2quat, quatrotate, round0
%
% -----------------------------------------------------------
%
function GI = getillumination( alpha, beta, gamma)

% normal of flat and unrotated solar cell pointing to zenith
n = [0,0,1];

% define euler matrix for rotation 'ZYX'
eul = deg2rad([alpha -beta gamma]);
rotm = eul2quat(eul);

% rotate normal about alpha and beta
coord_rot = quatrotate(rotm,n);
n_new = [coord_rot(1), coord_rot(2), coord_rot(3)];

% define local coodinates
phi = linspace(-pi,pi,361);
theta = linspace(-pi/2,pi/2,181);

% make meshgrid with phi and theta
[P,T] = meshgrid(phi,theta);
R = ones(size(P));

% transform the angle grid to cartesian coordinates
[px, py, pz] = sph2cart(P,T,R);
p(:,:,1) = px;
p(:,:,2) = py;
p(:,:,3) = pz;

% reshape the rotated normal to fit the size of the defined grid p
n_new = repmat(reshape(n_new,[1 1 3]),[size(p,1) size(p,2) 1]);

% calculate the angle between the rotated normal and the unrotated
% reference system p, both in cartesian coordinates
A=acosd(sum(n_new.*p,3));

% for angles smaler than 90° fill matrix GI with ones, otherwise GI is zero
% ATTENTION: round is needed! otherwise 90 is smaller than 90 =)
GI = double(round(A,2)<89); 

% take front side of cell only
GI(1:91,:) = 0;
GI = flipud(GI);

% imagesc(-180:180,0:180,GI)

end

