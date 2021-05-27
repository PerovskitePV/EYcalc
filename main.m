%% #########################
%  ### ENERGY YIELD MAIN ###
%  #########################

% ### PATH ###
addpath(genpath(pwd));  

% ### DATABASE ###
% Use database to store simulations and load already simulated data
StoreInDatabase = false;
PathOpticsResults = 'C:\Users\Raphael\Desktop\tmpoptics\';  % Please change!
PathEYResults = 'C:\Users\Raphael\Desktop\tmpEY\';          % Please change!

% ### REFRACTIVE INDEX DATA ###
% Load n,k data from databse
if ~exist('IndRefr','var')
    IndRefr = struct;
    [IndRefr.nkdata,IndRefr.names]=xlsread('_RefractiveIndexLib.xlsx');
end



%% #########################
%  ### IRRADIANCE MODULE ###
%  #########################
% Note: Takes up to 20 min, only needs to be executed once on each system! 

% ### LOCATION ###
% Select location for the Energy Yield calculations
CodeLocation  = '722020TYA';     % Code to be looked up from \Irradiance\Dataset_TMY3\User Manual TMY3.pdf
AliasLocation = 'Miami';         % To be specified


% ### LOAD / CALCULATE IRRDIANCE ###
% Define the folder name of choosen location
FolderNameIrradiance = strcat( pwd, ['\Irradiance\Spectra_',num2str(CodeLocation),'_',num2str(AliasLocation)] );

% Simulate the irradiance data
if exist(FolderNameIrradiance, 'dir') ~= 7
    Irradiance(CodeLocation, AliasLocation);
end
% Load the irradiance data, if it has been calculated already  
if ~exist('irradiance','var')
    irradiance = load(['Irradiance/Spectra_',num2str(CodeLocation),'_',num2str(AliasLocation),'/Irr_spectra_clouds.mat']);
    load(['Irradiance/Spectra_',num2str(CodeLocation),'_',num2str(AliasLocation),...
        '/TMY3_',num2str(CodeLocation),'_',num2str(AliasLocation),'.mat']);
    irradiance.Data_TMY3 = Data_TMY3; clear Data_TMY3;
end




%% #####################
%  ### OPTICS MODULE ###
%  #####################

% ### INPUT OPTICS MODULE ###
% Define the layer stack (names are from the refractive index database) with corresponding layer thicknesses in nm.
% The polarization can be either 'mixed', 'TE or 'TM'. For EY calculations, use 'mixed'.
% The wavelength range and angle resolution defines the resolution for the TMM simulations.
% The morphology is defined facing upwards (direction of incident light) and needs to be set for any incoherent layer.
% In order to use rear-textured c-Si one needs to texture the last incohrent layer in the stack (e.g. air).
% Incohrent layers can be defined manually and are automatically assumed, if the layer thickness exceeds 5µm
% A bifacial treatment of the stack involves a second simulation for light incident from the rear side
% The names of the absorbers can be specified as cell. An empty cell leads to auto-detection of the absorbers
Stack = {'Air','MgF2','Glass1.5','ITOfront','SnO2','Pero1.62','SpiroOMeTAD',...
    'ITOfront','aSi(n)','aSi(i)','cSi','aSi(i)','aSi(p)','ITOfront','Air'};
LayerThickness = [inf,100,1E4,100,10,450,20,25,5,5,250E3,5,5,100,inf];
Morphology = {'Flat','RandomUpright','RandomUpright'};
Polarization   = 'mixed';
lambdaTMM      = 300:5:1200;  
AngleResolution = 5;
IncoherentLayers = {'Air','Glass','EVA','Encapsulation','PDMS','cSi'};
bifacial = false;
Absorbers = {};

% ### CALL OPTICS MODULE ###
optics = OpticsModule(IndRefr, Stack, LayerThickness, AngleResolution, Morphology, bifacial,...
    Polarization, lambdaTMM, PathOpticsResults, StoreInDatabase, IncoherentLayers, Absorbers);     


% ### PLOT OPTICS ###
Aall = squeeze(sum(optics.Absorptance(:,1,:)));

figure; 
plot(lambdaTMM, optics.A );hold on
plot(lambdaTMM, optics.R );
plot(lambdaTMM, optics.T ); 
plot(lambdaTMM, Aall,'--');
plot(lambdaTMM, Aall + optics.R + optics.T, ':' ); hold off
xlabel('Wavelength (nm)'); ylabel('R, A, T')



%% ########################
%  ### ELECTRICS MODULE ###
%  ########################

% ### INPUT ELECTRICS MODULE ###
% Either define the electrical parameters of the multijunction solar cell
% and single cells or load experimental IV data.
% In case of defined electrical paramters, temperature effects are taken
% into account by temperature coefficients for Jsc and Voc. Moreover, the
% nominal temperature of the sub cells can be defined.
% In case of experimental data, these temperure effects are not taken into
% account.
electrics.configuration = '2T';     % 2T, 3T, 4T, 2T exp, 3T exp, 4T exp
electrics.shunt = 'with';               % with, without
electrics.RshTandem = 1000;             % shunt resistance of tandem device
electrics.RsTandem = 3;                 % serial resistance of tandem device
electrics.Rsh = [1300, 1000];           % shunt resistance of n-th cell
electrics.Rs = [2, 1];                  % serial resistance of n-th cell
electrics.CE = [1, 1];                  % collection efficiency of n-th cell
electrics.j0 = [2.7e-18, 1e-12];        % reverse-blocking current of n-th cell
electrics.n = [1.1, 1];                 % ideality factor of n-th cell
electrics.Temp = [25, 25];              % temperature of cells (can also be n vectors)
electrics.NOCT = [48, 48];              % nominal temperature of n-th cell, if a number, Temp is overwritten
electrics.tcJsc = [0.0002, 0.00032];    % temperature coefficient of Jsc in K^-1 of n-th cell
electrics.tcVoc = [-0.002, -0.0041];    % temperature coefficient of Voc in K^-1 of n-th cell

% Load experimental IV data for 1sun illumination.
% For "2T exp", only the tandem IV "electrics.IVtandem" can be defined.
% For "4T exp", the single cells (limited two 2 cells) need to be defined.
% The IV variable has to be a two collumn array, where the first collumn is
% the voltage, the second collumn the current.

% electrics.IVtop = load();
% electrics.IVtop = load();
% electrics.IVtandem = load();



%% ###########################
%  ### ENERGY YIELD MODULE ###
%  ###########################

% ### INPUT ENERGY YIELD MODULE ###
% Define rotation/orientation of solar cell by SolarCellRotationAngle (0=north, 90=east, 180=south, 270=west) and 
% SolarCellTiltAngle (rotation angle in degree about new y axis (for SolarCellRotationAngle = 0, SolarCellTiltAngle > 0 
% the cell tilts to the southern hemisphere). For more help: see EnergyYield() function
SolarCellRotationAngle = 180;
SolarCellTiltAngle = 20;

% Define if tracking of the solar cell should be enabled/disabled. The following options are available:
%  0:  disable tracking
%  1:  1-axis non-tilted east-west
%  2:  2-axis tracking
%  3:  1-axis latitude-titled zenith rotation
%  4:  1-axis latitude-tiled seesaw rotation; limitation: rotation needs to be = 180°
tracking = 0;           

% Take albedo into account. If enabled albedo is taken into account for tilted solar cells enhancing the current generation for front
% side illumination. If a bifacial solar cell is simulated, albedo enhances light intensity coming from the rear.
albedo = 0;

% Select a ground type for albedo simulations. Choose a filename from the Ecospec librabry.
% By default only black and white are available. You can get the Ecospec librabry here:
% https://speclib.jpl.nasa.gov/
% groundtype = 'artificialwhite';
groundtype = 'artificialblack';


% ### CALL ENERGY YIELD MODULE ###
if length(SolarCellTiltAngle)==1
    EY = EnergyYield(irradiance, optics, electrics, SolarCellRotationAngle, SolarCellTiltAngle, tracking, albedo, groundtype, PathEYResults, StoreInDatabase); 
else
    [EYaoi, TandemPowerTotal] = sweepEY(irradiance, optics, electrics, SolarCellRotationAngle, SolarCellTiltAngle, tracking, albedo, groundtype, PathEYResults, StoreInDatabase);
end


% ### PLOT EY ###
% figure; 
% subplot(3,2,1); plot(EY.Power_Tandem);    hold on; plot(EY.Power);        hold off; ylabel('Power (W/m²)');     xlim([0 8760]);
% subplot(3,2,2); plot(100*EY.FF_Tandem);   hold on; plot(100*EY.FF);       hold off; ylabel('FF (%)');           xlim([0 8760]);
% subplot(3,2,3); plot(EY.Voc_Tandem);      hold on; plot(EY.Voc);          hold off; ylabel('Voc (V)');          xlim([0 8760]);
% subplot(3,2,4); plot(EY.Jsc);                                                       ylabel('Jsc (mA/cm²)');     xlim([0 8760]);
% subplot(3,2,5); plot(EY.VMPP_Tandem);     hold on; plot(EY.VMPP);         hold off; ylabel('VMPP (V)');         xlim([0 8760]);
% subplot(3,2,6); plot(EY.TempModule(:,1)); hold on; plot(EY.TempAmbient);  hold off; ylabel('Temperature (°C)'); xlim([0 8760]);

% figure;
% plot(SolarCellTiltAngle,TandemPowerTotal);
% xlabel('Tilt angle (°)')
% ylabel('Anual Energy Yield (kWhm^{-2}a^{-1})')




%% ###########################
%  ### DATA POSTPROCESSING ###
%  ###########################
 
% ...
