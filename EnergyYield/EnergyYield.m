% ------------- DOCUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function calculates the energy yield of a
%                         tandem solar cell. It combines the irradiance
%                         data calculated by the irradiance module and the
%                         optics derived by the optics module. This
%                         function will also use the electrical module to
%                         calculate electric properties of the considered
%                         device. Additionally it is possible to rotate the
%                         solar cell about two angles (alpha and beta).
%          ^ N            The initial configuration is shown on the left.
%          | x            alpha rotated left-handed about the normal of the
%     y    |              solar cell (z axis), which follows the sun.
% W <------|------> O     About the rotated y axis beta rotates right-handed.
%          |              alpha=180, beta=30 means therefore tilt to southern
%          |              hemisphere.
%          v S            alpha=90, beta=20 equals a tilted solar cell about
%                         20° which is facing eastern hemisphere.
%
%
% #INPUT:                 irradiance (struct) - incl. all irradiance data
%                         optics (stuct) - incl. all optics data
%                         electrics (stuct) - incl. all electrics parameters
%                         lambda (array)
%                         alpha (scalar)
%                         beta (scalar)
%                         timetoanalyze (string or array)
%
% #OUTPUT:                EY (struct)
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: trimirradiance, getillumination, rotatesunangle,
%                         predef, calctandemelectrics, 
%                         calcsingleelectrics
%
% -----------------------------------------------------------
%
function EY = EnergyYield(irradiance, optics, electrics, SolarCellRotationAngle, SolarCellTiltAngle, tracking, albedo, groundtype, PathEYResults, StoreInDatabase)

% Begin optics code
warning('off','backtrace'); % hide line info in warnings
disp('Running EnergyYield Module...');
tic


% Check if variables are defined properly
sz0 = logical([1,1,1,0,1,0,1,0,1,1,0,1,1]);
fn = fieldnames(orderfields(electrics));
for i=1:numel(fn)
    sz(i)=size(electrics.(fn{i}),2);
end
if ~contains(electrics.configuration,'exp') && ~all(sz(sz0) == sz(1))
    warning('Input arguements in electrics struct have not the same size! Please define electrical properties for all absorbers!');
    toc
    return;
elseif size(optics.AbsorberIndex,2) ~= sz(1)
    warning('Number of absorbers in optics struct does not match the number of absorbers defined by the electrics!');
    toc
    return;
end

if numel(SolarCellTiltAngle)>1
    warning('Number of tilt angles is >1. You probably want to sweep. Please use sweepEY() instead of EnergyYield() for this!');
    toc
    return;
end  
    
if strcmp(electrics.configuration,'2T exp') && ~isfield(electrics,'IVtandem')
    warning('Attention: You selected experimental 2T electrics. In order to use an experimental IV curve for the 2T, you need to provide "electrics.IVtandem".')
    return;
end

if strcmp(electrics.configuration,'4T exp') && ( ~isfield(electrics,'IVtop') || ~isfield(electrics,'IVbot') )
    warning('Attention: You selected experimental 4T electrics. In order to use experimental IV curves for the 4T, you need to provide "electrics.IVtop" and "electrics.IVbot".')
    return;
end


%% Check for database to load energy yield data
if StoreInDatabase
    % Generate a string of the defined energy yield input parameters
    EYstring = strjoin( {   strjoin( strsplit( strjoin( cellfun(@num2str, struct2cell(electrics), 'un',0) ),' '), ',' ),...  % electrics
                            char(optics.hash), ...   % optics as hash
                            strjoin( {num2str(SolarCellTiltAngle), num2str(tracking), num2str(albedo), groundtype}, ',' ) } );   % rest of EY       

    % Generate SHA1 hash of the defined parameters and set some filenames
    Hash = makehash(EYstring);
    StrLine = strcat(strjoin(strcat([Hash,', ']),''), EYstring);
    
    % Build file names for optics data & database
    DataEYResults = strcat(PathEYResults,Hash ,'.mat');
    DBEYResults = strcat(PathEYResults,'_EYdatabase.txt');

    % Check if data is already in database & if the file exists
    % if true, load the already calculated data and return
    if ~isfile(DBEYResults)
        fclose(fopen(DBEYResults, 'w'));
    end
    if contains(fileread(DBEYResults), StrLine) && exist(DataEYResults, 'file') == 2
        load(DataEYResults);
        disp('EY Calculations has been performed before. Data is loaded from database.');
        toc
        return;
    end
end


%% Define some physical constants in Si units
q = 1.6021766208e-19;   % C 
h = 6.62607004e-34;     % Js = m^2 kg / s
heV = 4.135667662E-15;  % eV*s
c = 299792458;          % m/s
const = q/(h*c)*1e-10;  % precalculate constant for later Jsc calculations
lambda = optics.lambda(1):optics.lambda(end);   % wavelength with 1nm spacing for EY code


%% Load Optics
% assign absorption of absorbers in the stack
% A{front,back}, if monofacial A = [#absorbers x 1], if bifacial A = [#absobers x 2]
A = cell(length(optics.AbsorberIndex),isfield(optics, 'Absorptance_back')+1);
A(:,1) = mat2cell(optics.Absorptance(optics.AbsorberIndex,:,:), ones(length(optics.AbsorberIndex),1) );

% check if bifacial data is available
if isfield(optics, 'Absorptance_back')  
    A(:,2) = mat2cell(optics.Absorptance_back(optics.AbsorberIndex_back,:,:), ones(length(optics.AbsorberIndex_back),1) );
end

% reshape matrices in A from [1 x AOI x lambdaTMM] to [AOI x lambdaTMM]
A = cellfun(@squeeze, A, 'uniform', 0); 

% interpolate data to 1 nm resolution
A = cellfun(@(x) interp1(optics.lambda,x',lambda,'linear','extrap')', A, 'uniform', 0);


%% Load Irradiance
% Assign unrotated azimuth and elevation of sun
thetasun0 = irradiance.Data_TMY3(:,7);
phisun0 = irradiance.Data_TMY3(:,8);

% Shift angles by 0.5h to account for the coarse binning
% this resolved the steep drop of EY in the evening for tracking
thetasun0 = (thetasun0 + circshift(thetasun0,1))/2;
phisun0 = (phisun0 + circshift(phisun0,1))/2;

% Assign the wavelengths from irradiance data 
w = irradiance.Irr_spectra_clouds_wavelength(1,:);

% Assign and define complete (full spectral resolution) irradiance
% right now the irradiance module will provide horizontal irradiance only,
% therefore we need to calculate the normal irradiance by ourselfs
% Be careful with the difference for diffuse and direct:
% For the direct part this can be done by dividing through the cosine of
% the suns' zenith angle. In case of Lambertian distributed diffuse
% irradiance the normal part is equal to direct*1/pi.
IrradianceDifH = irradiance.Irr_spectra_clouds_diffuse_horizontal;
IrradianceDifN = IrradianceDifH / pi;
IrradianceDirH = irradiance.Irr_spectra_clouds_direct_horizontal;
IrradianceDirN = (ceil(thetasun0)<90) .* IrradianceDirH ./ repmat(abs(cosd(thetasun0)),[1, length(w)]);

% Ambient temperature
TempAmbient = irradiance.Data_TMY3(:,14);

% To minimize the amount of memory needed clean irradiance structure
clear irradiance;

% Trim irradiance data to wavelength range simulated by the optics module
% the irradiance is given as direct and diffuse horizontal. Both variables
% are matrices with size: (hour)x(lambda).
IdifN = trimirradiance(lambda, IrradianceDifN, w);
IdifH = trimirradiance(lambda, IrradianceDifH, w);
IdirN = trimirradiance(lambda, IrradianceDirN, w);
IdirH = trimirradiance(lambda, IrradianceDirH, w);

% Load Albedo-Data for ground type & define the diffuse irradiance
% corresponding to the albedo reflectance, for which a Lambertian
% distribution will be assumed. IalbN is calculated analoge to IrradianceDirN
if size(A,2)>1 || albedo == 1   % if bifacial or albedo
    albedoR = loadalbedo(lambda, groundtype);
    IalbN = albedoR .* (IdifH + IdirH)/pi;
end

% To estimate the cells efficiency we calculate the global irradiance:
% Iglob = diffuse horizontal + direct normal. This data will not be rotated
% Attention: the full spectral range will be used, not the trimmed data!
Iglob = sum( IrradianceDifH + IrradianceDirN ,2);


%% Calculate the energy yield
% First we define a matrix GI for the local (and rotated) cell coordiante
% system. This matrix is one for theta and phi (local mathematical spherical
% coordinates) from wich light can reach the cells top surface and zero for
% combinations of those angles from where no light can shine on the cell.
% Later can apprear, when the sun hits the back due to some rotation or theta
% and phi values which would describe illumination from the ground (albedo).
% Direct illumination from the backside and albedo reflections from the
% backside and front side for bifacial modules is treated by GI_inv.
GI = getillumination( SolarCellRotationAngle, SolarCellTiltAngle, 0 );
GI_inv = 1 - GI;
    
% Next, we rotate the suns' theta and phi angles. This gives us the angles
% of the sun in our new rotated local solar cell coordinate system.
[phisun, thetasun] = rotatesunangle( SolarCellRotationAngle, SolarCellTiltAngle, 0, phisun0, thetasun0 );

% To calculate the Jsc we need to perform an intergration over the half
% sphere and over the wavelength. Therefore we need the differentials and
% the integration variables.
theta = deg2rad(linspace(0,89,90));
dtheta = deg2rad(1);
dphi = deg2rad(1);
dlambda = lambda(2)-lambda(1);


% This function predefins the EY structure with some result variables
EY  = predef(size(A));

% predefine some tmp variables
alpha = zeros(8760,1);
beta = alpha;
gamma = alpha;
S = zeros(8760,1);

%% Loop over all hours of the considered year
for j=1:8760
    
    % If tracking is enabled GI will be re-calculated each j
    if tracking     % don't trust the algorithms without a second check 
        switch tracking         
            case 1  % 1-axis non-tilted east-west
                alpha(j) = 90;
                xsun = sind(phisun0(j))*cosd(90-thetasun0(j)); % west east cmponent
                zsun = sind(90-thetasun0(j)); % up-down component
                beta(j) = atand(xsun/zsun);
                gamma(j) = 0;
                
            case 2  % 2-axis tracking
                alpha(j) = phisun0(j);
                beta(j) = thetasun0(j);
                gamma(j) = 0;
                
            case 3  % 1-axis latitude-titled zenith rotation
                alpha(j) = phisun0(j);
                beta(j) = SolarCellTiltAngle;
                gamma(j) = 0;
                
            case 4  % 1-axis latitude-tiled seesaw rotation
                [ alpha(j), beta(j), gamma(j) ] = seesawtilt( SolarCellTiltAngle, SolarCellRotationAngle, phisun0(j), thetasun0(j) );
                
            case 5  % intensity tracking (not fully debuged yet) 
                sc=abs(sum(IdirN(j,:),2)-sum(IdifN(j,:),2))./(sum(IdirN(j,:),2)+sum(IdifN(j,:),2));
                alpha(j) = phisun0(j);
                beta(j) = sc*thetasun0(j);
        end
        if thetasun0(j)>=90
            beta(j) = 0;
            gamma(j) = 0;
        end
        
        [phisun(j), thetasun(j)] = rotatesunangle( alpha(j), beta(j), gamma(j), phisun0(j), thetasun0(j) );
        GI = getillumination( alpha(j), beta(j), gamma(j) );
        GI_inv = 1 - GI; 
        
    end
    
    % To model the module temperature we need to quantify the incident
    % insolation S
    S(j) = sum(IdirN(j,:) .* cosd(thetasun(j)) + IdifN(j,:) * sum( sum(ones(90,361),2)' .* sin(theta) .* cos(theta) ) * dphi * dtheta ) * dlambda;

    
    % Calculate the suns' phi and theta angle to be used as an index in GI
    idx_phisun = wrapTo360(round(phisun(j),0)+1);
    idx_thetasun = (round(thetasun(j),0)+1);
    
    % todo: check JscDirect < 1
    
    for k = 1:size(A,1)     % loop over number of absorber layers
        % direct light:
        if idx_thetasun<=90 && GI(idx_thetasun, idx_phisun) == 1	% direct hitting front for monofacial
            direct = (A{k,1}(idx_thetasun,:).*IdirN(j,:) * (electrics.CE(k) * lambda)')*cosd(thetasun(j));
        elseif size(A,2)>1 && 181-idx_thetasun<90 && thetasun0(j)<=90	% direct hitting back of bifacial
            direct = (A{k,2}(181-idx_thetasun,:).*IdirN(j,:)*(electrics.CE(k) .* lambda)')*cosd(181-thetasun(j));
        else	% direct not hitting anywhere
            direct = 0;
        end
        EY.JscDirect(j,k) = dlambda *  const * direct;
        
        % diffuse light:
        % monofacial w/o albedo
        tmp_diffuse = ( IdifN(j,:).* lambda )' * ( sum(GI(1:90,:),2)' .* sin(theta) .*cos(theta) );
        diffuse = sum( sum( electrics.CE(k) * (A{k,1}' .* tmp_diffuse) ));
        EY.JscDiffuse(j,k) = dphi * dtheta * dlambda * const * diffuse;
        
        if size(A,2)>1  % bifacial w/o albedo
            tmp_diffuse_back  = ( IdifN(j,:).* lambda )' * ( sum(GI_inv(1:90,:),2)' .* sin(theta) .*cos(theta) );
            diffuse = diffuse + sum( sum( electrics.CE(k) * (A{k,2}' .* tmp_diffuse_back) ));
            EY.JscDiffuse(j,k) = dphi * dtheta * dlambda * const * diffuse;
        end
        
        if albedo  % if albedo is enabled
            % monofacial
            tmp_alb  = ( IalbN(j,:).* lambda )' * ( sum(GI_inv(1:90,:),2)' .* sin(theta) .*cos(theta) );
            alb = sum( sum( electrics.CE(k) * (A{k,1}' .* tmp_alb) ));
            EY.JscAlbedo(j,k) = dphi * dtheta * dlambda * const * alb;
            
            if size(A,2)>1   % bifacial
                tmp_alb_back  = ( IalbN(j,:).* lambda )' * ( sum(GI(1:90,:),2)' .* sin(theta) .*cos(theta) );
                alb = alb + sum( sum( electrics.CE(k) * (A{k,2}' .* tmp_alb_back) ));
                EY.JscAlbedo(j,k) = dphi * dtheta * dlambda * const * alb;
            end
        end
    end
    EY.Jsc(j,:) =  EY.JscDirect(j,:) + EY.JscDiffuse(j,:) + EY.JscAlbedo(j,:);

end    % end of loop over time


%% Calculation of the module temperature based on the NOCT 
% (compare "Applied Photovoltaics" 3rd edition, page 79
if isnumeric(electrics.NOCT)
    electrics.Temp = zeros(8760,size(A,1));
    if length(electrics.NOCT)~=size(A,1)    % if same NOCT for all absobers
        electrics.NOCT = electrics.NOCT * ones(1,size(A,1));
    end
    for k=1:size(A,1)
        electrics.Temp(:,k) = TempAmbient + ( electrics.NOCT(k) - 20 ) / 800 .* S;
    end
end
EY.TempAmbient = TempAmbient;
EY.TempModule = electrics.Temp;

%% Calculate the electric characteristic of the solar cell(s) / tandem
% This depends on the electrics.configuration
if strcmp(electrics.configuration,'3T')
    warning('Although the 3T analytical model was selected, the 4T electrics model is used as a first approximation.');
end
switch electrics.configuration
    case '2T'
        % call tandem electrics with current matching
        [EY.Voc_Tandem, EY.FF_Tandem, EY.Power_Tandem, EY.JMPP_Tandem, EY.VMPP_Tandem] = calctandemelectrics(electrics, min(EY.Jsc,[],2));
                
        for k=1:size(A,1)
            % calculate single cells in tandem device in 2T
            [EY.Voc(:,k), EY.FF(:,k), EY.Power(:,k), EY.JMPP(:,k), EY.VMPP(:,k)] = ...
                calcsingleelectrics( min(EY.Jsc,[],2), electrics.Rs(k), electrics.Rsh(k), electrics.j0(k), electrics.n(k), electrics.tcJsc(k), electrics.tcVoc(k), electrics.Temp(:,k), electrics.shunt );

            % calculate single cells separate
            [EY.Voc_SJ(:,k), EY.FF_SJ(:,k), EY.Power_SJ(:,k), EY.JMPP_SJ(:,k), EY.VMPP_SJ(:,k)] = ...
                calcsingleelectrics( EY.Jsc(:,k), electrics.Rs(k), electrics.Rsh(k), electrics.j0(k), electrics.n(k), electrics.tcJsc(k), electrics.tcVoc(k), electrics.Temp(:,k), electrics.shunt );
        end
        
%   case '3T'
        % .. to be done.
        % .. at the moment the 4T model is used! This is indeed not correct, but is a solid approxmimation
        
    case {'4T','3T'}
        for k=1:size(A,1)
            % calculate single cells sperate
            [EY.Voc_SJ(:,k), EY.FF_SJ(:,k), EY.Power_SJ(:,k), EY.JMPP_SJ(:,k), EY.VMPP_SJ(:,k)] = ...
                calcsingleelectrics( EY.Jsc(:,k), electrics.Rs(k), electrics.Rsh(k), electrics.j0(k), electrics.n(k), electrics.tcJsc(k), electrics.tcVoc(k), electrics.Temp(:,k), electrics.shunt );
        end
        
        % calculate the power conversion efficiency of the tandem in 4T
        EY.Power_Tandem = sum(EY.Power_SJ,2);
    
    % Use experimental IV data. Attention, no temperature scaling is used.
    case '2T exp'  % 2T experimental IV data is limited to single tandem IV
        % smooth IV
        V = electrics.IVtandem(:,1);
        I = smooth(electrics.IVtandem(:,2));
        
        % format IV data to be in 1st quadrant
        [V,newIdx] = sort(V);
        I = I(newIdx);
        if I(1)<0; I = -I; end
            
        % sweep over hours
        for k=1:length(EY.Jsc)
            % scale current by Jsc
            I = I - ( I(find(diff(V>=0),1+1))  - min(EY.Jsc(k,:)) );   % I - ( I(V=0) - Jsc )
            P = (I .* V) * 10;
            
            [EY.Power_Tandem(k),~] = max(P,[],1);
        end

    case '4T exp'  % 4T experimental IV data is limited to 4T cell with 2 single cells
        % smooth IV
        Vtop = electrics.IVtop(:,1);
        Vbot = electrics.IVbot(:,1);
        Itop = smooth(electrics.IVtop(:,2));
        Ibot = smooth(electrics.IVbot(:,2));
        
        % format IV data to be in 1st quadrant
        [Vtop,newIdx] = sort(Vtop);
        Itop = Itop(newIdx);
        [Vbot,newIdx] = sort(Vbot);
        Ibot = Ibot(newIdx);
        
        if Itop(1)<0; Itop = -Itop; end
        if Ibot(1)<0; Ibot = -Ibot; end
        
        % sweep over hours
        for k=1:length(EY.Jsc)
            % scale current by Jsc
            Itop = Itop - ( Itop(find(diff(Vtop>=0),1+1))  - EY.Jsc(k,1) );   % I - ( I(V=0) - Jsc )
            Ibot = Ibot - ( Ibot(find(diff(Vbot>=0),1+1))  - EY.Jsc(k,2) );   % I - ( I(V=0) - Jsc )
            
            Ptop = (Itop .* Vtop) * 10;
            Pbot = (Ibot .* Vbot) * 10;
            
            [EY.Power_SJ(k,1),~] = max(Ptop,[],1);
            [EY.Power_SJ(k,2),~] = max(Pbot,[],1);
            EY.Power_Tandem = sum(EY.Power_SJ,2);
        end        

end


            
%% Define some other usefull variables
EY.TandemPCE = EY.Power_Tandem ./ Iglob;
EY.TandemPCE(isnan(EY.TandemPCE)) = 0;

EY.PCE = EY.Power ./ Iglob;
EY.PCE(isnan(EY.PCE)) = 0;

% add irradiance to EY structure
EY.IrradianceDifH = sum(IrradianceDifH,2);
EY.IrradianceDifN = sum(IrradianceDifN,2);
EY.IrradianceDirH = sum(IrradianceDirH,2);
EY.IrradianceDirN = sum(IrradianceDirN,2);

% add absorption Atop and Abot to EY
EY.A = A;

% do some post processing of hourly resolved data
EY.TandemPowerTotal = sum(EY.Power_Tandem)/1000;
EY.TandemPCEmean = mean(EY.TandemPCE);

% Albedo yes/no
EY.albedo = albedo;



%% Save data
if StoreInDatabase
    % write hash to database
    fileID = fopen(DBEYResults,'a');
    fprintf(fileID, strcat(StrLine,'\n'));
    fclose(fileID);
    
    % Save the optics data to a file
    save(DataEYResults,'EY');
end


disp(char(strcat('...finished after',{' '},num2str(toc),{' '},'seconds')));

end







