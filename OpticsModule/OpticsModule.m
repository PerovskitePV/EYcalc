% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function calculates the optics of an arbitrary layer stack consisting of optically thick and thin layers,
%                         whereas the optically thick layers are treated incohrent and can be textured. 
%                         The coherent treatment of light propagation is calculated with the transfer matrix method (TMM). For the TMM
%                         calculations, we adapted parts of the TMM code of Steven Byrnes.
%                         The optically thick layers are treated by simple calculation of the absorption in the dielectric layers via 
%                         the Beer-Lambert law. Those optically thick layers can additionally be textured. These textures can be handled 
%                         using geometrical ray-tracing as described by Baker-Finch and McIntosh. For this, path data is extracted from 
%                         OPAL2. Own path data can be used as well!
%                           1. The full stack is devided in incoherent and coherent partial layer stacks
%                           2. Cohrent stack (with incohrent boundaries are calculated with TMM
%                           3. Incoherent layers are treated with Lambert Beer and their interfaces with Fresnel
%                           4. PathRAT function redistributes the R,T,A matrices [aoi × lambda] to tensors [aoi × exit angle × lambda]
%                           5. R,T,A properties of Partial stacks are interativly connected via AddIncoherent or AddCoherentIncoherent
%
% #INPUT:                 IndRefr           - Refractive index structure with two fields: nkdata, names
%                         Stack             - Cell with layer names 
%                         LayerThicknesses  - Array with layer thicknesses in nm
%                         AngleResolution   - Angle resolution of TMM calculations (should not be >10°)
%                         Morphology        - Angle resolution of TMM calculations (should not be >10°)
%                         bifacial          - Boolean for bifacial option
%                         pol               - Polarization of optical simulations
%                         lambda            - Wavelength array in nm
%                         PathOpticsResults - Path to store and check for optics results
%                         StoreInDatabase   - Boolean for storing simulations in a database
%                         IncoherentLayers  - Cell with strings layer names, which should be treated incoherently - can be empty cell
%                         Absorbers         - Cell with names of the absorbers in the stack - can be empty cell
%                         
%
% #OUTPUT:               optics struct:
%                          Reflectance: [aoi × exit angle × lambda]
%                          Absorptance: [layer × aoi × lambda]
%                          Transmittance: [aoi × lambda]
%                          lambda: [1 × lambda]
%                          AbsorberIndex: [array]
%                          hash: [string]
% %
% #SAVED DATA:           The database file: "_OPTICSdatabase.txt" and corresponding *.mat files of the simulation results
% 
% #REQUIRED SUBFUNCTIONS: AddIncoherent, AddIncoherentCoherent, coh_tmm, coh_tmm_reverse, flipmorph, Fresnel, interpTMM, Lambert, 
%                         LoadRefrIndex, makehash, PathRAT 
%
% #ADD COMMENTS:         Angle redistribution with PathRAT does lead to underestimated A in textured interfaces close to the bandgap
%                        e.g. absorption in single or double-sided textured Si is slightly off. Easy fix: increase Si thickness by ~ *2
%                        in case no reflecting layer (e.g. metal) is on rear texture. See documentation for details.
%                        pathRATnew() tries to improve this by calculating the transmittance angles and following the bounces.
%                        This however, this only slightly improves the issue and slows down the simulation time.  
%                     
% -----------------------------------------------------------
%
function optics = OpticsModule(IndRefr,Stack,LayerThicknesses,AngleResolution,Morphology,bifacial,pol,lambda,PathOpticsResults,StoreInDatabase,IncoherentLayers,Absorbers)

% Begin optics code
warning('off','backtrace'); % hide line info in warnings
disp('Running Optics Module...');
tic

%% Check input for problems in the input to this function
if length(LayerThicknesses)~=numel(Stack)
    error('foo:bar',['Number of Layers (',num2str(numel(Stack)),...
        ') does not match the number of thicknesses (',num2str(length(LayerThicknesses)),')!' ]);
    return; %#ok<UNRCH>
end
if ~exist('Absorbers','var')
    if ~exist('IncoherentLayers','var')
        warning('No incoherent layers are defined! The default will be used: air, glass, eva, encapsulation and pdms');
        IncoherentLayers = {'Air','Glass','EVA','Encapsulation','PDMS'};    % define default incohrent layers
    end
    warning('Absorbers in current stack are not defined. No worries, we try to take care of that.');
    Absorbers = {};
elseif isempty(Absorbers)
    warning('Absorbers in current stack are not defined. No worries, we try to take care of that.');
    Absorbers = {};
end
    

%% Make hash(es) and get/make filenames and filename of database
[DataOpticsResults,DBOpticsResults,StrLine,Hash]=gennames;


%% Check for database to load optics data
if StoreInDatabase
    % Create database file, if not there
    if ~isfile(DBOpticsResults)
        fclose(fopen(DBOpticsResults, 'w'));
    end
    
    if contains(fileread(DBOpticsResults), StrLine(1,:)) && exist(DataOpticsResults(1,:), 'file') == 2
        load(DataOpticsResults(1,:));
        disp('Simulation was performed before. Data is loaded from database.');
    else
        optics = runoptics;
        optics.hash = Hash(1,:);
        saveresults(DataOpticsResults(1,:),StrLine(1,:))
    end
    
    if bifacial == 1
        if contains(fileread(DBOpticsResults), StrLine(2,:)) && exist(DataOpticsResults(2,:), 'file') == 2
            load(DataOpticsResults(2,:));
            disp('Simulation was performed before. Data is loaded from database.');
        else
            Stack = flip(Stack);
            LayerThicknesses = flip(LayerThicknesses);
            Morphology = flip(Morphology);
            optics = runoptics;
            optics.hash = Hash(2,:);
            saveresults(DataOpticsResults(2,:),StrLine(2,:))
        end
        optics_back = optics; 
        load(DataOpticsResults(1,:));
        optics.hash = optics_back.hash+optics.hash;
        optics.Reflectance_back = optics_back.Reflectance;
        optics.Absorptance_back = optics_back.Absorptance;
        optics.Transmittance_back = optics_back.Transmittance;
        optics.AbsorberIndex_back = numel(Stack) - optics.AbsorberIndex - 1;
    end
else    
    optics = runoptics;
    optics.hash = Hash(1,:);
    if bifacial == 1
        Stack = flip(Stack);    
        LayerThicknesses = flip(LayerThicknesses);
        Morphology = flip(Morphology);
        
        optics_back = runoptics;
        
        optics.hash = Hash(1,:)+Hash(2,:);
        optics.Reflectance_back = optics_back.Reflectance;
        optics.Absorptance_back = optics_back.Absorptance;
        optics.Transmittance_back = optics_back.Transmittance;
        optics.AbsorberIndex_back = numel(Stack) - optics.AbsorberIndex - 1;
    end
end

disp(char(strcat('...finished after',{' '},num2str(toc),{' '},'seconds')));




%% Core optics function
function optics = runoptics

    %% Define some stuff
    MorphologyBack = flipmorph(Morphology);     % flip morphologies for backward simulations
    TMMAngleRange=[0:AngleResolution:88,89];    % define angle range for TMM
    AngleRange=0:89;                            % define angle range for outpt
    ThicknessThreshold = 5E3;                   % threshold of cohrent layer in nm


    %% Predefine variables and load material data
    n = zeros(size(Stack,2),size(lambda,2));
    for index = 1:size(Stack,2)
        n(index,:) = LoadRefrIndex(Stack{index},lambda,IndRefr.nkdata,IndRefr.names);
    end
    

    %% Analyze the stack
    % Check for a match in Stack
    Match = or( cell2mat(cellfun(@(x) contains(x, IncoherentLayers), Stack, 'UniformOutput', 0)), LayerThicknesses>ThicknessThreshold );

    % Define temporary variables to track the positions of incohrent or
    % coherent layers as well as to organize the stack splitting
    %   pos_i       positions of incohrent layers
    %   pos_ii      position of incohrent/incohrent interface in Rf,Tf,..
    %   pos_c       position of cohrent layers
    %   pairs_i     position of two successive incohrent layers
    %   pairs_tmm   position of incohrent layers forming TMM boundaries
    %   pairs_all   position of each partial stack
    pos_i = 1:length(Match);
    pos_i = pos_i(Match(1:end));
    pairs_i = [strfind(Match, [1,1]);strfind(Match, [1,1])+1];
    pos_ii = 1:length(pos_i); 
    pos_ii = pos_ii(((pos_i(2:end))-pos_i(1:end-1))==1);    pos_ii(isempty(pos_ii))=0; % set 0, if not found
    pos_c = 1:length(pos_i);
    pos_c = pos_c(logical(1-sum((pos_c-pos_ii')==0,1)));
    pos_c = pos_c(1:end-1); 
    pairs_tmm = reshape(sort([strfind(Match, [1 0]), strfind(Match, [0 1])+1]),2,[]);
    pairs_all = sort([pairs_tmm,pairs_i],2);

    % Check if each incohrent layer has a morphology
    if numel(Morphology) ~= length(pos_i)-1
        error('foo:bar',['The number of assigned morphologies (',num2str(numel(Morphology)), ...
            ') does not match the number of incoherent layers (',num2str(length(pos_i)-1),')!\n\n', ...
            'Note: For each incoherent layer (except the first) a textured needs to be defined.']);
        return; %#ok<UNRCH>
    end

    % Preallocate the cells to store "RAT" data of TMM or Fresnel inside
    [Rf, Af, Tf, Rb, Ab, Tb] = deal(cell(size(pos_i,2)-1,1));  % initialize R,T fw and bw


    %% Calculate Lambert Beer absorption in incoherent layers
    Aincoherent = cell(length(pos_i)-2,1); % initialize cell
    for i=2:length(pos_i)-1 % for each (i)<last
        Aincoherent{i-1} = Lambert(n(pos_i(i),:),LayerThicknesses(pos_i(i)),lambda);
        %     figure; imagesc(Aincoherent{i});
    end


    %% Calculate Fresnel reflection and transmittance at the interface of two incoherent layers
    % Fresnel(n1,n2,pol,lambda) takes care of absorbing media by switching to TMM in case imag(n1)~=0
    % This is essential for e.g. air/cSi - here Fresnel() would fail in the short wvl region
    % due to strong absorption coefficient in cSi.
    for i=1:size(pairs_i,2)
        j=pos_ii(i);    % index of position where two incohrent layers are
        [Rf{j},Tf{j}] = Fresnel((n(pairs_i(1,i),:)),(n(pairs_i(2,i),:)),pol,lambda);   % R,T forward
        [Rb{j},Tb{j}] = Fresnel((n(pairs_i(2,i),:)),(n(pairs_i(1,i),:)),pol,lambda);   % R,T backward
        [Af{j},Ab{j}] = deal(zeros([1 size(Rf{j})]));
        %     figure; imagesc(Rf{j});
        %     figure; imagesc(Tf{j});
    end


    %% Calculate TMM for TMM stacks
    % Preallocate cells
    [Rf_TMM_TE, Af_TMM_TE, Tf_TMM_TE, Rb_TMM_TE, Ab_TMM_TE, Tb_TMM_TE] = deal(cell(size(pairs_tmm,2),1));
    [Rf_TMM_TM, Af_TMM_TM, Tf_TMM_TM, Rb_TMM_TM, Ab_TMM_TM, Tb_TMM_TM] = deal(cell(size(pairs_tmm,2),1));

    % Loop over number of TMM stacks
    for i=1:size(pairs_tmm,2)
        j=pos_c(i);    % index of position where two TMM layer stack is
        b = pairs_tmm(1,i):pairs_tmm(2,i);  % Layer range for TMM stack
        c = LayerThicknesses(b);            % Thicknesses of i-th TMM stack
        c([1,end]) = Inf;                   % Change incoherent boundary to Inf

        if strcmp(pol,'TE') || strcmp(pol,'TM')        % Check pol and calculate TMM
            for aoi=1:length(TMMAngleRange)     % Loop over angle of incidence
                [Rf{j}(aoi,:),Af{j}(:,aoi,:),Tf{j}(aoi,:),~] = coh_tmm(pol, n(b,:).', c, TMMAngleRange(aoi), lambda);
                [Rb{j}(aoi,:),Ab{j}(:,aoi,:),Tb{j}(aoi,:),~] = coh_tmm_reverse(pol, n(b,:).', c, TMMAngleRange(aoi), lambda);
            end

        elseif strcmp(pol,'mixed')
            for aoi=1:length(TMMAngleRange)
                % TMM for TE - fw and bw
                [Rf_TMM_TE{i}(aoi,:),Af_TMM_TE{i}(:,aoi,:),Tf_TMM_TE{i}(aoi,:),~] = coh_tmm('TE', n(b,:).', c, TMMAngleRange(aoi), lambda);
                [Rb_TMM_TE{i}(aoi,:),Ab_TMM_TE{i}(:,aoi,:),Tb_TMM_TE{i}(aoi,:),~] = coh_tmm_reverse('TE', n(b,:).', c, TMMAngleRange(aoi), lambda);
                % TMM for TM - fw and bw
                [Rf_TMM_TM{i}(aoi,:),Af_TMM_TM{i}(:,aoi,:),Tf_TMM_TM{i}(aoi,:),~] = coh_tmm('TM', n(b,:).', c, TMMAngleRange(aoi), lambda);
                [Rb_TMM_TM{i}(aoi,:),Ab_TMM_TM{i}(:,aoi,:),Tb_TMM_TM{i}(aoi,:),~] = coh_tmm_reverse('TM', n(b,:).', c, TMMAngleRange(aoi), lambda);
            end
            % Average TE and TM polarization in order to obtain mixed polarization
            Rf{j} = ( Rf_TMM_TE{i} + Rf_TMM_TM{i} ) /2;
            Af{j} = ( Af_TMM_TE{i} + Af_TMM_TM{i} ) /2;
            Tf{j} = ( Tf_TMM_TE{i} + Tf_TMM_TM{i} ) /2;
            Rb{j} = ( Rb_TMM_TE{i} + Rb_TMM_TM{i} ) /2;
            Ab{j} = ( Ab_TMM_TE{i} + Ab_TMM_TM{i} ) /2;
            Tb{j} = ( Tb_TMM_TE{i} + Tb_TMM_TM{i} ) /2;

            % Interpolate the data to 1 deg of incidence
            [Rf{j}, Af{j}, Tf{j}, ~] = interpTMM(Rf{j}, Af{j}, Tf{j}, n(pairs_tmm(1,i),:), n(pairs_tmm(2,i),:), c, TMMAngleRange, lambda);
            [Rb{j}, Ab{j}, Tb{j}, ~] = interpTMM(Rb{j}, Ab{j}, Tb{j}, n(pairs_tmm(1,i),:), n(pairs_tmm(2,i),:), c, TMMAngleRange, lambda);
        end
    end
    % for i=1:size(Rf,1)
    %     figure; imagesc(Tf{i})
    % end

    
    %% Calculate RAT of textured interfaces
    % There are two possibilities. PathRAT() and PathRATnew(). The first
    % one, is a rough and stable function. Since, the path data do not
    % provide transmittance angles, this version produces underestimation
    % of absorption e.g. in textured c-Si near the bandgap.
    % In order to avoid this, the angles are explicitly calculated.
    % PathRATnew() tries to improve this. However, altough the angular
    % distribution matches almost the raytracing results, still there is an
    % offset. This function, still needs to be improved.
    for i=1:length(pos_i)-1
        if ischar(Morphology{i}) && ~strcmp(Morphology{i},'4n^2')
            
            % Calculate forward RAT
            [Rf{i},Af{i},Tf{i}] = PathRAT(Morphology{i},Rf{i},Af{i},Tf{i},AngleRange,lambda,n(pairs_all(1,i),:),n(pairs_all(2,i),:));
%             [Rf{i},Af{i},Tf{i}] = PathRATnew(Morphology{i},Rf{i},Af{i},Tf{i},AngleRange,lambda,n(pairs_all(1,i),:),n(pairs_all(2,i),:));
            
            % Calculate backward RAT
            if i<length(pos_i)-1
                [Rb{i},Ab{i},Tb{i}] = PathRAT(MorphologyBack{i},Rb{i},Ab{i},Tb{i},AngleRange,lambda,n(pairs_all(2,i),:),n(pairs_all(1,i),:));
%                 [Rb{i},Ab{i},Tb{i}] = PathRATnew(MorphologyBack{i},Rb{i},Ab{i},Tb{i},AngleRange,lambda,n(pairs_all(2,i),:),n(pairs_all(1,i),:));
            
            % set backward RAT for last interface to zero
            else
                [Rb{i},Tb{i}] = deal(zeros(size(Rf{i})));
                Ab{i} = zeros(size(Af{i}));
            end
            
        end
    end

    
    %% Connect the partial stacks
    for i=1:length(pos_i)-2 % loop over all partial stacks, maybe replace by: size(pairs_all,2)-1 
        % break >last< loop if silicon is calculated with a light trapping model
        if i==length(pos_i)-2 && (~ischar(Morphology{end}) || strcmp(Morphology{end},'4n^2'))
            break;
        end
        % x + i ...+i
        if numel(Match(pos_i(i):pos_i(i+2))) == 3 || sum(Af{i+1}(:))==0
            if i==1 % if this is the first occurance
                [Rftot,Aftot,Tftot,Rbtot,Abtot,Tbtot] = AddIncoherent(Rf{i},Af{i},Tf{i},Rb{i},Ab{i},Tb{i},Aincoherent{i},Rf{i+1},Rb{i+1},Tf{i+1},Tb{i+1});
            else
                [Rftot,Aftot,Tftot,Rbtot,Abtot,Tbtot] = AddIncoherent(Rftot,Aftot,Tftot,Rbtot,Abtot,Tbtot,Aincoherent{i},Rf{i+1},Rb{i+1},Tf{i+1},Tb{i+1});
            end
        % x + i+c ...+i
        else
            if i==1 % if this is the first occurance  
                [Rftot,Aftot,Tftot,Rbtot,Abtot,Tbtot] = AddIncoherentCoherent(Rf{i},Af{i},Tf{i},Rb{i},Ab{i},Tb{i},Aincoherent{i},Rf{i+1},Af{i+1},Tf{i+1},Rb{i+1},Ab{i+1},Tb{i+1});
            else
                [Rftot,Aftot,Tftot,Rbtot,Abtot,Tbtot] = AddIncoherentCoherent(Rftot,Aftot,Tftot,Rbtot,Abtot,Tbtot,Aincoherent{i},Rf{i+1},Af{i+1},Tf{i+1},Rb{i+1},Ab{i+1},Tb{i+1});
            end
        end
        % reshape transmittance data in last loop
        if i == (length(pos_i)-2)
            Tftot = reshape(sum(Tftot,2),size(Tftot,1),size(Tftot,3));
        end
    end

    % if only coherent stack with incoherent bondary is given
    if size(Rf,1)==1
        Rftot = Rf{1};
        Aftot = Af{1};
        Tftot = reshape(sum(Tf{1},2),size(Tf{1},1),size(Tf{1},3));
    end

    
    % Calculate light trapping model for Silicon:
    % if only Si (e.g.)
    if ( strcmp(Morphology{end},'4n^2') || (isnumeric(Morphology{end}) && ~strcmp(Morphology{end},'4n^2')) ) && ~exist('Rftot','var')
        Rftot = Rf{1};
        Aftot = Af{1};
        Tftot = Tf{1};
    end
    if strcmp(Morphology{end},'4n^2')
        Z = 4*real(n(pairs_all(1,end),:)).^2;
        Aftot(end+1,:,:) = reshape(sum(Tftot,2),size(Tftot,1),size(Tftot,3)) .* repmat(1-exp(-(4*pi./lambda)...
                .*imag(n(pairs_all(1,end),:)).*LayerThicknesses(pairs_all(1,end)).*Z),[90 1]);
        Tftot = reshape(sum(Tftot,2),size(Tftot,1),size(Tftot,3)) - reshape(Aftot(end,:,:),size(Aftot,2),size(Aftot,3));    
    elseif isnumeric(Morphology{end}) && ~strcmp(Morphology{end},'4n^2')
        if ~isnumeric(Morphology{end})
            Morphology{end} = str2double(Morphology{end});
        end
        Z = Morphology{end}*ones(1,length(lambda));
        Aftot(end+1,:,:) = reshape(sum(Tftot,2),size(Tftot,1),size(Tftot,3)) .* repmat(1-exp(-(4*pi./lambda)...
                .*imag(n(pairs_all(1,end),:)) .* LayerThicknesses(pairs_all(1,end)).*Z),[90 1]);
        Tftot = reshape(sum(Tftot,2),size(Tftot,1),size(Tftot,3)) - reshape(Aftot(end,:,:),size(Aftot,2),size(Aftot,3));
    end




    %% Put everything into the optics struct
    % Results are converted to single precission - this is sufficient!
    % This saves 50% of storage and speeds up saving/loading.
    optics.Reflectance = single(Rftot);
    optics.Absorptance = single(Aftot);
    optics.Transmittance = single(Tftot);
    optics.lambda = lambda;

    
    %% auto detect absorber layers, in case nothing is defined
    if isempty(Absorbers)
        A = sum(sum(Aftot(:,1,:),3),2); A=A/max(A);
        if length(A)<3
            [~,optics.AbsorberIndex]=max(A);
        else
            [~,optics.AbsorberIndex]=findpeaks([A;0]','Threshold',1/3);
        end
    else       
        optics.AbsorberIndex = find(matches(Stack,Absorbers) == 1);
        optics.AbsorberIndex = optics.AbsorberIndex - (length(Stack) - size(Aftot,1) - 1); % correct for offset     
    end

    
    %% Save vectors for normal incidence and absorbers only
    optics.A = squeeze(optics.Absorptance(optics.AbsorberIndex,1,:))';
    optics.R = squeeze(sum(optics.Reflectance(1,:,:),2));
    optics.T = squeeze(optics.Transmittance(1,:))';
    
    
end



%% #########################
% Helper functions below
%# #########################

%% Generate filename function
function [DataOpticsResults,DBOpticsResults,StrLine,Hash] = gennames
    % Generate a string of the defined optics input parameters
    opticstring = strjoin({strjoin(Stack),...
        strjoin(strsplit(num2str(LayerThicknesses)),' '),...
        strjoin(Morphology),...
        pol,...
        num2str(lambda(1)), num2str(lambda(2)-lambda(1)), num2str(lambda(end)),...
        num2str(AngleResolution),...
        strjoin(Absorbers),...
        strjoin(IncoherentLayers)},', ');

    % Generate SHA1 hash of the defined parameters and set some filenames
    Hash(1,:) = makehash(opticstring);
    StrLine = strcat(strjoin(strcat([Hash(1,:),', ']),''), opticstring);

    % Build file names for optics data & database
    DataOpticsResults = strcat(PathOpticsResults,Hash(1,:) ,'.mat');
    DBOpticsResults = strcat(PathOpticsResults,'_OPTICSdatabase.txt');

    if bifacial == 1
        % Generate a string of the defined optics input parameters
        opticstring = strjoin({strjoin(flip(Stack)),...
            strjoin(strsplit(num2str(flip(LayerThicknesses))),' '),...
            strjoin(flip(Morphology)),...
            pol,...
            num2str(lambda(1)), num2str(lambda(2)-lambda(1)), num2str(lambda(end)),...
            num2str(AngleResolution),...
            strjoin(Absorbers),...
            strjoin(IncoherentLayers)},', ');

        % Generate SHA1 hash of the defined parameters and set some filenames
        Hash(2,:) = makehash(opticstring);
        StrLine(2,:) = strcat(strjoin(strcat([Hash(2,:),', ']),''), opticstring);

        % Build file names for optics data & database
        DataOpticsResults(2,:) = strcat(PathOpticsResults,Hash(2,:) ,'.mat');
    end
end
   

%% Save data
function saveresults(fname,strname)
    % write hash to database
    fileID = fopen(DBOpticsResults,'a');
    fprintf(fileID, strcat(strname,'\n'));
    fclose(fileID);

    % Save the optics data to a file
    save(fname,'optics');
end




end

