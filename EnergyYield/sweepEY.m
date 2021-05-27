% #DESCRIPTION:           This function sweeps the EY simulation for various SolarCellTiltAngle
%
% #INPUT:                 same as EnergyYield()
%                         SolarCellTiltAngle can be a vector instead of a scalar
%
% #OUTPUT:                EYaoi - struct containing EY data for all SolarCellTiltAngle
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% -----------------------------------------------------------

function [EYaoi, TandemPowerTotal] = sweepEY(irradiance, optics, electrics, SolarCellRotationAngle, SolarCellTiltAngle, tracking, albedo, groundtype, PathEYResults, StoreInDatabase)
    
    EYaoi = cell(length(SolarCellTiltAngle),1);
    TandemPowerTotal = zeros(length(SolarCellTiltAngle),1);
    
    parfor i=1:length(SolarCellTiltAngle)
        EYaoi{i} = EnergyYield(irradiance, optics, electrics, SolarCellRotationAngle, SolarCellTiltAngle(i), tracking, albedo, groundtype, PathEYResults, StoreInDatabase);
        EYaoi{i}.SolarCellTiltAngle = SolarCellTiltAngle(i);
    end

    for i=1:length(SolarCellTiltAngle)
        TandemPowerTotal(i) = EYaoi{i}.TandemPowerTotal;
    end

end