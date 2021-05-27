% ------------- DOCUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:           Calculates RAT along the paths of a textured 
%                         effective interface (e.g. air/ARC/glass) for each angle of incidence
% #INPUT:                 Morphology: Surface morphology of the Si solar cell (string)
%                         InterfaceReflection: Reflection of flat effective interface
%                         InterfaceAbsorption: Absorption of the flat effective interface
%                         InterfaceTransmission: Transmission of the flat effective interface
%                         AngleRange: Angle of incidence (vector)
%                         thickness: Layer thicknesses of effective interface (array)
%                         lambda: Wavelengths (vector)
%                         n1,n2: Refractive indices of starting and end layer (array)
% #OUTPUT:                TensorReflection: Total reflection of the effective interface (tensor)
%                         TensorAbsorption: Total absorption of the effective interface(tensor)
%                         TensorTransmission: Total transmission of the effective interface (tensor)
% #SAVED DATA:            -
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          Available surface morphologies: 
%                           Flat
%                           RegularUpright
%                           RandomUpright
%                           Inverted
%                         
% -----------------------------------------------------------

function [TensorReflection,TensorAbsorption,TensorTransmission] = PathRAT(Morphology,InterfaceReflection,InterfaceAbsorption,InterfaceTransmission,AngleRange,lambda,n1,n2)
     
    % Preallocate path data
    pathFractions = cell(1,90);
    anglesIntersect = cell(1,90);
    exitAngle = cell(1,90);
    firstFacet = cell(1,90);
    facetNormals = cell(1,90);
    
    num_layers = size(InterfaceAbsorption,1);
    num_angles = length(AngleRange);
    num_lambdas = length(lambda);
        
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
    
    % Calculate RAT textured surface
    TensorReflection=zeros(num_angles,num_angles,num_lambdas);
    TensorAbsorption=zeros(num_layers,num_angles,num_lambdas);
    TensorTransmission=zeros(num_angles,num_angles,num_lambdas);
    for alphaInIdx = 1:num_angles
        PathReflection=zeros(length(pathFractions{alphaInIdx}),num_lambdas);
        PathAbsorption=zeros(num_layers,length(pathFractions{alphaInIdx}),num_lambdas);
        PathTransmission=zeros(length(pathFractions{alphaInIdx}),num_angles,num_lambdas);
        for i = 1:length(pathFractions{alphaInIdx}) % Paths
            alphaOutIdx=exitAngle{alphaInIdx}(i)+1;
            v = [sind(alphaInIdx-1) 0 -cosd(alphaInIdx-1)];
            n = facetNormals{firstFacet{alphaInIdx}(i)};
            vtrans = snell3D(v,n,real(n1),real(n2));
            alphaTransIdx = floor(round(abs(atand(vtrans(1,:)./vtrans(3,:))),10))+1;    % round needed because of flow(x.99999) = x instead of x+1!
            if max(alphaTransIdx > 90) == 1
                alphaTransIdx(alphaTransIdx > 90 & alphaTransIdx <= 180) = alphaTransIdx(alphaTransIdx > 90 & alphaTransIdx <= 180) - 90; % Check for total internal reflection. This should be handled correctly by TMM.
                alphaTransIdx(alphaTransIdx > 180 & alphaTransIdx <= 270) = 270 - alphaTransIdx(alphaTransIdx > 180 & alphaTransIdx <= 270); % Check for total internal reflection. This should be handled correctly by TMM.
                alphaTransIdx(alphaTransIdx > 270 & alphaTransIdx <= 360) = alphaTransIdx(alphaTransIdx > 180 & alphaTransIdx <= 270) - 270; % Check for total internal reflection. This should be handled correctly by TMM.
            end
            for j = 1:length(anglesIntersect{alphaInIdx}(1,:)) % Bounces per path
                if j==1
                    PathReflection(i,:) = InterfaceReflection(anglesIntersect{alphaInIdx}(i,j)+1,:);
                    PathAbsorption(:,i,:) = InterfaceAbsorption(:,anglesIntersect{alphaInIdx}(i,j)+1,:);
                    % Assign first transmission of first order to known
                    % transmittance angle
					for k = 1:num_lambdas
                        PathTransmission(i,alphaTransIdx(k),k) = InterfaceTransmission(anglesIntersect{alphaInIdx}(i,j)+1,k);
                    end
                elseif j>1
                    if isnan(anglesIntersect{alphaInIdx}(i,j))
                        break;
                    end
                    PathAbsorption(:,i,:) = PathAbsorption(:,i,:)+reshape(repmat(PathReflection(i,:),[num_layers 1]),num_layers,1,num_lambdas)...
                        .*InterfaceAbsorption(:,anglesIntersect{alphaInIdx}(i,j)+1,:);
                    PathTransmission(i,:,:) = PathTransmission(i,:,:)...
                         + reshape(repmat(PathReflection(i,:).*InterfaceTransmission(anglesIntersect{alphaInIdx}(i,j)+1,:)/num_angles, num_angles,1),1,num_angles,num_lambdas);
%                         + reshape( (sind((0:89)))' .* ( PathReflection(i,:).*InterfaceTransmission(anglesIntersect{alphaInIdx}(i,j)+1,:)/sum(sind((0:89))) ),1,num_angles,num_lambdas); % https://www.comsol.com/release/5.3/ray-optics-module
                     PathReflection(i,:) = PathReflection(i,:).*reshape(InterfaceReflection(anglesIntersect{alphaInIdx}(i,j)+1,:),1,num_lambdas);
                end
            end
            TensorTransmission(alphaInIdx,:,:) = TensorTransmission(alphaInIdx,:,:)...
                + pathFractions{alphaInIdx}(i)*PathTransmission(i,:,:);
            TensorReflection(alphaInIdx,alphaOutIdx,:) = TensorReflection(alphaInIdx,alphaOutIdx,:)...
                + reshape(pathFractions{alphaInIdx}(i)*PathReflection(i,:),1,1,num_lambdas);
        end
        
       % Calculate total A of top cell for each wavelength
       TensorAbsorption(:,alphaInIdx,:) = sum(repmat(pathFractions{alphaInIdx}.',[num_layers 1 num_lambdas]).*PathAbsorption,2);
       
    end
end