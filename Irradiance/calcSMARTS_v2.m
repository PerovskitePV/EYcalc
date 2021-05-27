% ------------- DOKUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:          Executes the SMARTS Code to calculate the clear sky irradiance
%
% #INPUT:                Code_location,Alias_location
%
% #OUTPUT:               /
%
% #SAVED DATA:           Saves clear Sky Irradiance 'Irr_spectra_clear_sky.mat'.
%                        FORMAT: Saved as *.mat file containing three matrixes (365x24, length([280:1:4000)).  
%                        1st matrix: 'Irr_spectra_clear_sky_wavelength' 
%                        2nd matrix: 'Irr_spectra_clear_sky_direct_normal'
%                        3rd matrix: 'Irr_spectra_clear_sky_diffuse_horizontal'
%                        4th matrix: 'Irr_spectra_clear_sky_direct_horizontal'
%
% #REQUIRED SUBFUNCTIONS: /
%
% #ADD COMMENTS:   The SMARTS code is a software acessible via NREL. Here
%                  we provide a convenient way to execute it.
%                  
%                  Attention: This function writes a lot of files to the
%                  disk! Depending on your system this will be
%                  significantly quicker than calcSMARTS(), which
%                  overwrites the input file for the SMARTS calculation for
%                  each hour of the year. Here, all input files are
%                  generated first (~4000-5000 txt files).
%
%-----------------------------------------------------------
%
function calcSMARTS_v2(Code_location,Alias_location)

% Define variables
NAME_location = [Code_location,'_',Alias_location];
folder_target_data = ['Irradiance\Spectra_',NAME_location];
load([folder_target_data,'\TMY3_',NAME_location,'.mat']);
location_SmartS = 'Irradiance\Code_SMARTS_295_PC';

% Predefine vectors and matrixes for data
Irr_spectra_clear_sky_wavelength         = 280:4000;  % Wavelength in nm
Irr_spectra_clear_sky_direct_normal      = zeros(365*24,length(280:4000),'single');
Irr_spectra_clear_sky_diffuse_horizontal = zeros(365*24,length(280:4000),'single');
Irr_spectra_clear_sky_direct_horizontal  = zeros(365*24,length(280:4000),'single');


location.zenith   = Data_TMY3(:,7);
location.azimuth  = Data_TMY3(:,8);
location.altitude = Data_TMY3(:,22);
location.longitude= Data_TMY3(:,20);
location.latitude = Data_TMY3(:,21);

specs.Irr_ext_hor_surf    = Data_TMY3(:,9);   % Hourly extraterrestrial horizontal surface  Wh/m2
specs.Irr_ext_norm_sun    = Data_TMY3(:,10);  % Hourly extraterrestrial radiation normal to the sun  Wh/m2
specs.Irr_hori_surface    = Data_TMY3(:,11);  % Global horizontal irradiance Wh/m2
specs.Irr_direct_normal   = Data_TMY3(:,12);  % Direct normal irradiance
specs.Irr_diff_horizontal = Data_TMY3(:,13);  % Diffuse horizontal irradiance Wh/m2
specs.dry_bulb_T          = Data_TMY3(:,14);  % Dry-bulb temperature  C
specs.humidity_relative   = Data_TMY3(:,15);  % relative humuduty  %
specs.pressure            = Data_TMY3(:,16);  % Station pressure mbar
specs.precitable_H2O      = Data_TMY3(:,17);  % preciptable water 1/cm
specs.AOD                 = Data_TMY3(:,18);  % aoreosol optical depth [unitless]   used for turbidity
specs.albedo              = Data_TMY3(:,19);  % albedo [unitless]

    
% Change to SMRATS folder
here = pwd;
cd([location_SmartS,'\']);

% Define imput for SMARTS
for k = 1:8760
    
    if specs.Irr_ext_norm_sun(k)~=0
        filename = [num2str(k,'%04d'),'.inp.txt'];

        Output.C1   =   ['''',num2str(k),''''];                          % Card 1:   Comment
        Output.C2   =   '1';                                                                                        % Card 2:    Controll - site pressure - details next 
        Output.C2a  =   [num2str(specs.pressure(k)),' ',num2str(location.altitude(k)/1000),' 0'];                          % Card 2a:   'pressure[mB] altitude[km] height[km]'
        Output.C3   =   '1';                                                                                        % Card 3:    IATMOS - 0: no reference available 1: reference 
        Output.C3a  =   '''USSA''';                                                                                  % Card 3a:   if 1: 'USS'; reference from 1976 
        Output.C4   =   '0';                                                                                        % Card 4:    Water vapour 
        Output.C4a  =   num2str(specs.precitable_H2O(k));                                                               % Card 4a:   0: precitable water in cm
        Output.C5   =   '1';                                                                                         % Card 5:    Ozone
        Output.C6   =   '1';                                                                                         % Card 6:    Gaseous adsorption
        Output.C7   =   '370';                                                                                       % Card 7:    Carbon dioxide 370
        Output.C7a  =   '0';                                                                                         % Card 7a:   extraterrestical spectrum ... 0 - Gueymard 2004
        Output.C8   =   '''S&F_RURAL''';                                                                             % Card 8:    Aerosol Model   'S&F_RURAL'  Rurual Shettle and Fenn
    %   Output.C8a  =   '0'                                                                                         % Card 8a:  
        Output.C9   =   '0';                                                                                         % Card 9:    Turbidity 0: AOD at 500nm / we use average value though
        Output.C9a  =   num2str(specs.AOD(k));                                                                          % Card 9a:   value [unitless]
        Output.C10  =   '-1';                                                                                        % Card 10:   Albedo -1 fixed value from 
        Output.C10a =   num2str(specs.albedo(k));                                                                       % Card 10a:  value [unitless]
        Output.C10b =   '0';                                                                                         % Card 10b:  no tilted albedo
        Output.C11  =   ['280 4000 1.0 ', num2str(specs.Irr_ext_norm_sun(k))];                                          % Card 11:   Min & max wavelengths; sun-earth distance correction; solar constant  
        Output.C12  =   '2';                                                                                         % Card 12:   Number of Ouputs 
        Output.C12a =   '280 4000 1';                                                                              % Card 12:   [extraterrestial, dir. norm. irr, dir. tilted irr, diff. tilted irr, global tilted irr., exp direct w. circumsolar, exp diff. irr 
        Output.C12b =   '10';                                                                                        % Card 12:   Number of Ouputs 
        Output.C12c =   '1 2 3 4 5 6 7 8 9 10';                                                                    % Card 12:   [extraterrestial irr, dir. norm. irr, dir. horizontal irr, diff. horzontal irr, global horizontal irr, direct horizontal irr 
        Output.C13  =   '1';                                                                                         % Card 13:   circumsolar - NOT RELEVANT 
        Output.C13a =   '0 5 0';                                                                                     % Card 13a:  circumsolar - NOT RELEVANT 
        Output.C14  =   '0';                                                                                         % Card 14:   smoothin filter - NOT RELEVANT 
        Output.C15  =   '0';                                                                                         % Card 15:   extra illuminance - NOT RELEVANT 
        Output.C16  =   '0';                                                                                         % Card 16:   extra UV calcs below 280 nm - NOT RELEVANT 
        Output.C17  =   '0';                                                                                         % Card 17:   Solar Geometry 

        if location.zenith < 90
            Output.C17a = [num2str(location.zenith(k)),' ',num2str(location.azimuth(k))];                                % Card 17:   Zenith - Azimu   
        else
            Output.C17a = [num2str(90-1E-3),' ',num2str(location.azimuth(k))];                                        % Card 17:   Zenith - Azimu   
        end
        Out =  {Output.C1;Output.C2;Output.C2a;Output.C3;Output.C3a;Output.C4;Output.C4a;Output.C5;Output.C6;Output.C7;...
            Output.C7a;Output.C8;Output.C9;Output.C9a;Output.C10;Output.C10a;Output.C10b;Output.C11;Output.C12;Output.C12a;...
            Output.C12b;Output.C12c;Output.C13;Output.C13a;Output.C14;Output.C15;Output.C16;Output.C17;Output.C17a};    

        writecell(Out,filename)
        
        parameters = {'n', '.', num2str(k,'%04d'), 'y'}';
        writecell(parameters,['para',num2str(k,'%04d'),'.txt']);
    end
  
end

% Call SMARTS
parfor k=1:8760
    if specs.Irr_ext_norm_sun(k)~=0
        system( sprintf('smarts295.exe < %s > NUL', ['para',num2str(k,'%04d'),'.txt']) );
        delete(['para',num2str(k,'%04d'),'.txt']);
        delete([num2str(k,'%04d'),'.inp.txt']);
        delete([num2str(k,'%04d'),'.out.txt'])
    end
end


% Load results
for k=1:8760
    if specs.Irr_ext_norm_sun(k)~=0      
        Data_import = importdata([num2str(k,'%04d'),'.ext.txt'],' ',1);
        
        delete([num2str(k,'%04d'),'.ext.txt']);
        
        Irr_spectra_clear_sky_direct_normal(k,:)      = interp1(Data_import.data(:,1),Data_import.data(:,3), [280:1:4000]');  % Direct_normal Irradiance in W/m2/nm
        Irr_spectra_clear_sky_diffuse_horizontal(k,:) = interp1(Data_import.data(:,1),Data_import.data(:,4), [280:1:4000]');  % Diffuse horizontal Irradiance in W/m2/nm
        Irr_spectra_clear_sky_direct_horizontal(k,:)  = interp1(Data_import.data(:,1),Data_import.data(:,6), [280:1:4000]');  % Direct horizontal Irradiance in W/m2/nm
    end
end

cd(here);



save([folder_target_data,'\Irr_spectra_clear_sky.mat'], 'Irr_spectra_clear_sky_wavelength', 'Irr_spectra_clear_sky_direct_normal', 'Irr_spectra_clear_sky_diffuse_horizontal', 'Irr_spectra_clear_sky_direct_horizontal');

% save([folder_target_data,'\Irr_spectra_clear_sky_wavelength'], 'Irr_spectra_clear_sky_wavelength')
% save([folder_target_data,'\Irr_spectra_clear_sky_direct_normal'], 'Irr_spectra_clear_sky_direct_normal')
% save([folder_target_data,'\Irr_spectra_clear_sky_diffuse_horizontal.mat'], 'Irr_spectra_clear_sky_diffuse_horizontal')
% save([folder_target_data,'\Irr_spectra_clear_sky_direct_horizontal.mat'], 'Irr_spectra_clear_sky_direct_horizontal')

fclose('all');


end
