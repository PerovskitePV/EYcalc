% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function loads and interpolates Albedo Data
%                         for the selected ground type. After interpolation
%                         of the loaded data, the resulting spectrum fits
%                         to the wavelength range of lambda defined in
%                         main.m. If necessary, the spectrum is
%                         extrapolated with a constant value to the
%                         required range of wavelength.
%
% #INPUT:                 lambda (vector)
%                         groundtype (string)
%
% #OUTPUT:                albedoR (vector)
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% -----------------------------------------------------------
%
function albedoR = loadalbedo(lambda, groundtype)

    fileID = fopen(['EnergyYield/LibraryEcospec/',groundtype,'.txt']);
    albedodata = sortrows(cell2mat(textscan(fileID,'%f %f','HeaderLines',21)));
    fclose(fileID);
    albedodata(:,1) = 1000 * albedodata(:,1);
    albedodata(:,2) = 1/100 * albedodata(:,2);
    albedoR = interp1(albedodata(:,1),albedodata(:,2),lambda);
    for o = 1:length(albedoR)
        if albedoR(o) < 0
            albedoR(o) = 0;
        end
        if isnan(albedoR(o)) && albedodata(o,1) < mean(lambda)
            albedoR(o) = albedodata(1,2);
        elseif isnan(albedoR(o)) && albedodata(o,1) > mean(lambda)
            albedoR(o) = albedodata(size(albedodata,1),2);
        end
    end

end