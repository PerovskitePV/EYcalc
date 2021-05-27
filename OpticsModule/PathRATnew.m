% ------------- DOCUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           Calculates RAT along the paths of a textured
%                         effective interface (e.g. air/ARC/glass) for each angle of incidence
% #INPUT:                 Morphology: Surface morphology of the Si solar cell (string)
%                         Ri: Reflection of flat effective interface
%                         Ai: Absorption of the flat effective interface
%                         Ti: Transmission of the flat effective interface
%                         AngleRange: Angle of incidence (vector)
%                         thickness: Layer thicknesses of effective interface (array)
%                         lambda: Wavelengths (vector)
%                         n1,n2: Refractive indices of starting and end layer (array)
% #OUTPUT:                R: Total reflection of the effective interface (tensor)
%                         A: Total absorption of the effective interface(tensor)
%                         T: Total transmission of the effective interface (tensor)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          Available surface morphologies:
%                           Flat
%                           RegularUpright
%                           RandomUpright
%                           Inverted
% -----------------------------------------------------------

function [R,A,T] = PathRATnew(Morphology,Ri,Ai,Ti,AngleRange,lambda,n1,n2)
% only works with mat files, where all facets are added

% todo:
%   code still is not much better than before..
%   check azimuth


% Preallocate path data
pathFractions = cell(1,90);         % path fractions - probability for each path
anglesIntersect = cell(1,90);       % intersection angles for each path and each aoi
exitAngle = cell(1,90);             % final exit angle of path
firstFacet = cell(1,90);            % all facets for each path and each aoi (in early version was only first..)
facetNormals = cell(1,90);          % facet normals

num_layers = size(Ai,1);            % number of layers
num_angles = length(AngleRange);    % number of angles
num_lambdas = length(lambda);       % number of wavelenghts

% Load path data from mat files
load(strcat('pathFractions',Morphology,'.mat'));
load(strcat('anglesIntersect',Morphology,'.mat'));
load(strcat('exitAngle',Morphology,'.mat'));
load(strcat('firstFacet',Morphology,'.mat'));
load(strcat('facetNormals',Morphology,'.mat'));

% Round path data
for idx=1:num_angles
    anglesIntersect{idx}=round(anglesIntersect{idx});
    anglesIntersect{idx}(anglesIntersect{idx}==90)=89;
    exitAngle{idx}=round(exitAngle{idx});
    exitAngle{idx}(exitAngle{idx}==90)=89;
end


% Preallocate R,A,T for the textured surface
R=zeros(num_angles,num_angles,num_lambdas);
A=zeros(num_layers,num_angles,num_lambdas);
T=zeros(num_angles,num_angles,num_lambdas);

% Calculate the TIR angle
TIR = ceil(real(asind(real(n2)./real(n1))));

% AOI
for alphaInIdx = 1:num_angles
    PathReflection=zeros(length(pathFractions{alphaInIdx}),num_angles,num_lambdas);     % path x aoi x lambda
    PathAbsorption=zeros(num_layers,length(pathFractions{alphaInIdx}),num_lambdas);     % layers x path x lambda
    PathTransmission=zeros(length(pathFractions{alphaInIdx}),num_angles,num_lambdas);   % path x aoi x lambda
    
    % Paths
    for pathIdx = 1:length(pathFractions{alphaInIdx})
        alphaOutIdx(1) = exitAngle{alphaInIdx}(pathIdx)+1;      % "normal" exit ray of this path
        v = [sind(alphaInIdx-1) 0 -cosd(alphaInIdx-1)];         % Incident ray vector
        n = facetNormals{firstFacet{alphaInIdx}(pathIdx)};      % First facet normal
        vtrans = snell3D(v,n,real(n1),real(n2));                % Transmitted vector
        alphaTransIdx = anglevec3d(vtrans, [0 0 -1], true)+1;   % Transmission angle (index)

        % Bounce
        for bounceIdx = 1:length(anglesIntersect{alphaInIdx}(1,:)) % Bounces per path
            
            if bounceIdx==1 % First bounce
                PathReflection(pathIdx,alphaOutIdx(1),:) = Ri(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,:);   % Reflection at #1 bounce
                PathAbsorption(:,pathIdx,:) = Ai(:,anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,:);              % Absorption at #1 bounce
                
                % Transmission
                for k = 1:num_lambdas   % Sweep over all wavelenghts
                    
                    % First bounce transmitts in case of n1 < n2
                    if n1(k)<n2(k) && isreal(vtrans(:,k))
                        PathTransmission(pathIdx,alphaTransIdx(k),k) = ...                      
                            Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);   % Transmission for intersection angle at transmission angle
                    
                    % First bounce transmitts in case of n1 > n2: intersection angle < TIR, z component of ray vector > z component of normal
                    elseif n1(k)>n2(k) && anglesIntersect{alphaInIdx}(pathIdx,bounceIdx) < TIR(k) && abs(v(3)) > abs(n(3)) && isreal(vtrans(:,k))
                        PathTransmission(pathIdx,alphaTransIdx(k),k) = ...                      
                            Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);   % Transmission for intersection angle at transmission angle

                    % First bounce transmitts in case of n1 > n2: intersection angle < TIR, z component of ray vector < z component of
                    % normal, angle between ray vector and opposite normal > 90°
                    elseif n1(k)>n2(k) && anglesIntersect{alphaInIdx}(pathIdx,bounceIdx) < TIR(k) && abs(v(3)) < abs(n(3)) && anglevec3d(-v,n.*[-1 -1 1],false)+1 > 90 && isreal(vtrans(:,k))
                        t = (vtrans(:,k)' - 2 * dot(n.*[-1 -1 1],vtrans(:,k)) * (n.*[-1 -1 1]))';       % Transmitted ray vector after reflection at opposite facet
                        
                        alphaExtraTIdx = anglevec3d(t,[0 0 -1],false)+1;                                % Transmission angle still inside n1
                        if alphaExtraTIdx > 50 || ( abs(t(3,:)) > n(:,3) && alphaExtraTIdx > 50 )
                            t = (t' - 2 * dot(n,t) * n)';
                            alphaExtraTIdx = anglevec3d(t,[0 0 -1],false)+1;                                % Transmission angle still inside n1
                        end
                        

                        PathTransmission(pathIdx,alphaExtraTIdx,k) = ...        
                            Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);                     % Add transmission at new angle
                        
                        % transmission into n2 -> contribution to refelctance?
                        % check!
                        % need to scale Ti.. -> ray split.
                        %
                        
                    else 
                        % Adjust transmission, in case above rules fail
                        if PathReflection(pathIdx,alphaOutIdx(1),k) ~= 1
                            PathTransmission(pathIdx,alphaTransIdx(k),k) = Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);
                        end
                        % no light is transmitted. Nothing to do here.
                    end
                end
                
                
            elseif bounceIdx>1 % other bounces
                if isnan(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx))
                    break;
                end
                
   %% old
%                 PathAbsorption(:,pathIdx,:) = PathAbsorption(:,pathIdx,:) + reshape(repmat(PathReflection(pathIdx,alphaOutIdx(1),:),[num_layers 1]),num_layers,1,num_lambdas)...
%                      .* Ai(:,anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,:); 
% 
%                 PathTransmission(pathIdx,:,:) = PathTransmission(pathIdx,:,:)...
%                     + reshape(repmat(sum(squeeze(PathReflection(pathIdx,:,:))).*Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,:)/num_angles, num_angles,1),1,num_angles,num_lambdas);
%                 
%                 PathReflection(pathIdx,:,:) = PathReflection(pathIdx,:,:).*repmat( reshape(Ri(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,:),1,1,num_lambdas), 1,num_angles);
                
   %% new             
    
                % Absorption: previous path absorption + ( "path refelction of this path with exit angle: alphaOutIdx" * "interface
                % absorption for intersection angle at current bounce" )
   % !!!        % alphaOutIdx might need to be modified if multiple reflection angles exits
                PathAbsorption(:,pathIdx,:) = PathAbsorption(:,pathIdx,:) + reshape(repmat(PathReflection(pathIdx,alphaOutIdx(end),:),[num_layers 1]),num_layers,1,num_lambdas)...
                    .* Ai(:,anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,:); 
                
                % next transmission angle
                v = v - 2 * dot(n,v) * n;                                               % Reflected light vector
                n = facetNormals{firstFacet{alphaInIdx}(pathIdx,bounceIdx)};            % Next/new facet normal
                vtrans = snell3D(v,n,real(n1),real(n2));                                % Transmitted ray vectors for each wavelengths of #j bounce!
                alphaTransIdx = anglevec3d(vtrans,[0 0 -1],true)+1;

                
                % sweep over wavelenghts
                for k = 1:num_lambdas
                    
                    if anglesIntersect{alphaInIdx}(pathIdx,bounceIdx) <= TIR(k) && vtrans(3,k) < 0 % ray points down

                       if alphaTransIdx(k) <= 70 && n1(k)<n2(k)
                           
                            PathTransmission(pathIdx,alphaTransIdx(k),k) = PathTransmission(pathIdx,alphaTransIdx(k),k)...
                                + PathReflection(pathIdx,alphaOutIdx(1),k) .* Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);
                            PathReflection(pathIdx,alphaOutIdx(1),k) = PathReflection(pathIdx,alphaOutIdx(1),k) .* Ri(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);

                        elseif  alphaTransIdx(k) > 70 && n1(k)>n2(k) % ray couples in again
                            
                           r = snell3D( vtrans(:,k)', n.*[-1 -1 1] ,real(n2(:,k)), real(n1(:,k)) );    % additional reflectance
                           alphaExtraRIdx1 = anglevec3d(real(r),n.*[-1 -1 1],false)+1;
                           
                           r = (r' - 2 * dot(n,r') * (n))';                                            % typicall if it coules in, the ray makes another TIR bounce
                           alphaExtraRIdx2 = anglevec3d(real(r),[0 0 1],false)+1;
                           
                           t = (vtrans(:,k)' - 2 * dot(n.*[-1 -1 1],vtrans(:,k)) * (n.*[-1 -1 1]))';   % transmittance
                           alphaExtraTIdx = anglevec3d(real(t),[0 0 -1],false)+1;
                           
                           PathTransmission(pathIdx,alphaExtraTIdx,k) = PathTransmission(pathIdx,alphaExtraTIdx,k)...
                                + PathReflection(pathIdx,alphaOutIdx(1),k) .* Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k) * Ri(alphaExtraRIdx1,k);
                            
                           PathReflection(pathIdx,alphaOutIdx(1),k) = PathReflection(pathIdx,alphaOutIdx(1),k) .* Ri(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k); 
                           PathReflection(pathIdx,alphaExtraRIdx2,k) = Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k) *  Ti(alphaExtraRIdx1,k);

                        end
                        
                    elseif ( anglesIntersect{alphaInIdx}(pathIdx,bounceIdx) <= TIR(k) && isreal(vtrans(3,k)) && vtrans(3,k) >= 0 ) || Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k) > 0 % ray points up (or horizontal)
                        
                        % first, adjust reflected ray
                        PathReflection(pathIdx,alphaOutIdx(1),k) = PathReflection(pathIdx,alphaOutIdx(1),k) * ...
                            Ri(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);
                        
                        % if ray is transmitted
                        if Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k) > 0
                            
                            r = snell3D( vtrans(:,k)', n.*[-1 -1 1] ,real(n2(:,k)), real(n1(:,k)) );    % additional reflectance
                            t = (vtrans(:,k)' - 2 * dot(n.*[-1 -1 1],vtrans(:,k)) * (n.*[-1 -1 1]))';   % transmittance
                            
                            alphaExtraRIdx = anglevec3d(real(r),[0 0 1],false)+1;
                            alphaExtraTIdx = anglevec3d(real(t),[0 0 -1],false)+1;
   
                            % ray gets bend upwards, so it will interact most likely again with the texture
                            % if so, calculate new angle and assume TIR
                            if alphaTransIdx(k) < alphaExtraTIdx
                                t = (t' - 2 * dot(n,t) * n)';
                                alphaExtraTIdx = anglevec3d(real(t),[0 0 -1],true)+1;
                            end
                            alphaIntersect = anglevec3d(n.*[-1 -1 1], real(vtrans(:,k)), true)+1;

                            if alphaExtraRIdx < 90  
                                PathReflection(pathIdx,alphaExtraRIdx,k) = PathReflection(pathIdx,alphaExtraRIdx,k) + ...
                                     Ti(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k) * Ti(alphaIntersect,k);
                                 %    (1 - PathReflection(pathIdx,alphaOutIdx(1),k)) * Ti(alphaIntersect,k);
                                
                                 
                                PathAbsorption(:,pathIdx,k) = PathAbsorption(:,pathIdx,k) + (1 - sum(PathReflection(pathIdx,:,k)) - sum(PathTransmission(pathIdx,:,k)) ) * Ai(:,alphaIntersect,k);
                               
                                PathTransmission(pathIdx,alphaExtraTIdx,k) = PathTransmission(pathIdx,alphaExtraTIdx,k) + ...
                                    (1 - sum(PathReflection(pathIdx,:,k)) - sum(PathTransmission(pathIdx,:,k)) - sum(PathAbsorption(:,pathIdx,k),1) );
                                   % (1 - sum(PathReflection(pathIdx,:,k)) - sum(PathTransmission(pathIdx,:,k)) ) * Ri(alphaIntersect,k);
                               % 1 ist falsch . fix vll so wie jetzt - aber dann fehlt immernoch die absorption.
                               
                                
                            else
                                PathTransmission(pathIdx,alphaExtraTIdx,k) = PathTransmission(pathIdx,alphaExtraTIdx,k) + ...
                                    (1 - sum(PathReflection(pathIdx,:,k)) - sum(PathTransmission(pathIdx,:,k)) );
                            end
                        end
                        
                        % debug
                        % figure(80); imagesc(squeeze(sum(PathTransmission,1)))
                        % figure(90); imagesc(squeeze(sum(PathReflection,1)))
                        % figure(70);plot(squeeze(sum(sum(PathAbsorption))));hold on;plot(squeeze(sum(PathTransmission(:,alphaInIdx,:))));plot(squeeze(sum(PathReflection(:,alphaInIdx,:))));hold off
           
                    else % no light is transmitted. Nothing to do here.
                        PathReflection(pathIdx,alphaOutIdx(1),k) = PathReflection(pathIdx,alphaOutIdx(1),k) * ...
                            Ri(anglesIntersect{alphaInIdx}(pathIdx,bounceIdx)+1,k);
                        
                        % A?
                        continue;
                    end
                    
                end
                
            end
        end
        

        if ~all( round(squeeze(sum(PathReflection(pathIdx,:,:),2)) + squeeze(sum(PathTransmission(pathIdx,:,:),2)) + squeeze(sum(PathAbsorption(:,pathIdx,:),1)),3) == 1)
            warning('R+T not one..')
        end
        
        % adjust transmittance and reflectance by path fraction
        T(alphaInIdx,:,:) = T(alphaInIdx,:,:) + pathFractions{alphaInIdx}(pathIdx) * PathTransmission(pathIdx,:,:);
        R(alphaInIdx,:,:) = R(alphaInIdx,:,:) + pathFractions{alphaInIdx}(pathIdx) * PathReflection(pathIdx,:,:);
        
    end
    
    % Calculate total A of top cell for each wavelength
    A(:,alphaInIdx,:) = sum(repmat(pathFractions{alphaInIdx}.',[num_layers 1 num_lambdas]).*PathAbsorption,2);
    
    if ~all( round( squeeze( sum(T(alphaInIdx,:,:),2) + sum(R(alphaInIdx,:,:),2) + sum(A(:,alphaInIdx,:),1) ),3) == 1)
        warning(['R+T+A at',num2str(alphaInIdx),' not one..'])
    end
        
end


end













