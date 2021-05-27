% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function calculates physical properties of
%                         tandem solar cell for a 2-terminal monolitic
%                         architecture based on electrical 1-diode-model
%                         without shunt resistance for each subcell. The
%                         influence of solar cell temperature is
%                         implemented using temperature coefficients for
%                         VOC and JSC.
%
% #INPUT:                 serial resistance: Rs
%                         reverse-blocking currents: j0_top, j0_bottom
%                         ideality factor: n_top, n_bottom
%                         short-circuit currents: jsc_top, jsc_bottom
%                         temperature coefficient of JSC: tcJSC_top,
%                         tcJSC_bottom
%                         temperature coefficient of VOC: tcVOC_top,
%                         tcVOC_bottom
%                         temperature of solar cell: Temp_top, Temp_bot
%
% #OUTPUT:                open circuit voltage: TandemVOC
%                         fill factor: TandemFillfactor
%                         electrical power: TandemP_el
%                         current at MPP: TandemJMPP
%                         voltage at MPP: TandemVMPP
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS: -
% -----------------------------------------------------------


function [TandemVOC,TandemFF,TandemP_el,TandemJMPP,TandemVMPP] = calctandemelectrics(electrics, jsc_RT)

% Number of subcells
k = max( [size(electrics.Rs,2), size(electrics.Rsh,2), size(electrics.j0,2), size(electrics.n,2)] );

% Thermal voltage at room temperature in V
Vth = 0.02569;

% Preallocate variables
[tmpJsc,tmpVoc,jsc,VOC_RT,VOC,thrs] = deal(zeros(size(jsc_RT,1),k));

for i=1:k
    % Temperature scaling factors
    tmpJsc(:,i) = (1 + electrics.tcJsc(i) * (electrics.Temp(:,i) - 25));
    tmpVoc(:,i) = (1 + electrics.tcVoc(i) * (electrics.Temp(:,i) - 25));
    
    % JSC for single cells
    jsc(:,i) = jsc_RT .* tmpJsc(:,i);
    
    % VOC for single cells
    thrs(:,i) = jsc(:,i) > 100*electrics.j0(i);

    % Calculate Voc at room temperature
    if strcmp(electrics.shunt,'without')
        VOC_RT(:,i) = electrics.n(i) * Vth * log(jsc(:,i)./electrics.j0(i));
    elseif strcmp(electrics.shunt,'with')
        VOC_RT(:,i) = jsc(:,i)/1000 * electrics.RshTandem - electrics.n(i) * Vth * ...
          lambertwlog(log( electrics.j0(i)/1000 * electrics.RshTandem) + electrics.RshTandem * (jsc(:,i) + electrics.j0(i))/(1000 * electrics.n(i) * Vth)  - log (electrics.n(i) * Vth) ) + ...
          electrics.j0(i)/1000 * electrics.RshTandem;
    end
    VOC(:,i) = thrs(:,i) .* VOC_RT(:,i) .* tmpVoc(:,i);
    VOC(isnan(VOC(:,i)),i) = 0;
end

% VOC tandem at chosen temperatures for top cell and bottom cell
TandemVOC = sum(VOC,2);

% Current density in mA/cm^2
jstart = -5;
jpoints = 401;
jstep = 1e-1;

j = (jstart:jstep:(jstart+jstep*jpoints));

V = zeros(size(jsc_RT,1),length(j),k);

for i=1:k
    if strcmp(electrics.shunt,'without')
        Vcond = ( repmat(jsc(:,i), [1 length(j)]) - repmat(j,[length(jsc(:,i)) 1]) )/electrics.j0(i);
        V(:,:,i) = electrics.n(i) * Vth * log((Vcond>0).*Vcond) - VOC_RT(:,i) + VOC(:,i);
        V(:,:,i) = max(V(:,:,i),0);
    elseif strcmp(electrics.shunt,'with')
        lambw = log( electrics.j0(i)/1000 * electrics.RshTandem ) + ( electrics.RshTandem * (-j + jsc(:,i) + electrics.j0(i))/(1000 * electrics.n(i) * Vth) ) - log ( electrics.n(i) * Vth );
        V(:,:,i) = - j/1000 * ( electrics.RsTandem + electrics.RshTandem) + jsc(:,i)/1000 * electrics.RshTandem - electrics.n(i) * Vth .* ...
            lambertwlog( lambw ) + electrics.j0(i)/1000 * electrics.RshTandem - VOC_RT(:,i) + VOC(:,i);
    end
end

% Tandem Voltage
V_tandem = (prod(V,3)~=0) .* ( sum(V,3) - repmat(j,[size(jsc,1) 1])/1000 * electrics.RsTandem );
V_tandem(isnan(V_tandem))=0;
V_tandem(isinf(V_tandem))=0;

% Tandem Power
Ptandem = j .* V_tandem;
[P_max,I] = max(Ptandem,[],2);

% Generate output
TandemJMPP = j(I)' .* prod(thrs,2);%[j] = mA/cm^2
for z=1:length(I)
    TandemVMPP(z,:) = V_tandem(z,I(z));%[V] = V
end
TandemJSC = min(jsc,[],2);%[j] = mA/cm^2
TandemFF = TandemJMPP .* TandemVMPP ./ (TandemJSC .* TandemVOC);%FF dimensionless
TandemFF(isnan(TandemFF))=0;
TandemP_el = P_max.*10;%[P] = W/m^2

end



