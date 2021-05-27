% ------------- DOKUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:          irradiance.m organises the derivation of the spectra
%                              1. Extract climate and meteo data from TMY3 dataset
%                              2. Calculate clear sky irradiace w. SMARTS
%                              3. Calcualte real spectra w. cloud model
%
% #INPUT:                INPUT_code_location, INPUT_ALIAS_location
%
% #OUTPUT:               /
%
% #SAVED DATA:           Subfunctions saves spectra in the folder ['Spectra_'Code_location'_'Alias_Location']
%                        Detailed description of data found in respective subfunctions
%
% #REQUIRED SUBFUNCTIONS: 
%         1. extractTMY3(Code_location,Alias_location)
%         2. calcSMARTS(Code_location,Alias_location)
%         3. simpleclouds(Code_location,Alias_location)
%
% #ADD COMMENTS:         
%   1. You might need to be admin to excecute the command line the calculate_spectra_from_SMARTS_code
%   2. Depending on the speed of your system, the pause time in "calculate_spectra_from_SMARTS_code" might need to be extended
%
%-----------------------------------------------------------
%
function Irradiance(Code_location,Alias_location)

% Extract the data from the TMY3 datasets
extractTMY3(Code_location,Alias_location)

% Calculate the irradiance spectra with help of the SMARTS code
calcSMARTS(Code_location,Alias_location)

% Enhance irradiance by a simply clouds model
simpleclouds(Code_location,Alias_location)

end
