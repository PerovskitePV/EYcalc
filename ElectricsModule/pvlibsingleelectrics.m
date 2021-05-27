% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%
% #DESCRIPTION:           This function calculates physical properties of one
%                         subcell based on electrical 1-diode-model without
%                         shunt resistance. The influence of solar cell
%                         temperature is implemented using temperature
%                         coefficients for VOC and JSC.
%
% #INPUT:                 serial resistance: Rs
%                         reverse-blocking current: j0
%                         ideality factor: n
%                         short circuit current: jsc
%                         temperature coefficient of JSC: tcJSC
%                         temperature coefficient of VOC: tcVOC
%                         temperature of solar cell: Temp
%
% #OUTPUT:                open circuit voltage: SingleVOC
%                         fill factor: SingleFillfactor
%                         electrical power: SingleP_el
%                         current at MPP: SingleJMPP
%                         voltage at MPP: SingleVMPP
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS: Test function using pvlib route to calculate IV curve
%                Attention: compared to pvlib, the temperature effects on
%                           Jsc and Voc are taken into account here.
% -----------------------------------------------------------

function [VOC,FF,P_el,JMPP,VMPP,j,V] = pvlibsingleelectrics(jsc_RT, Rs, Rsh, j0, n, tcJSC, tcVOC, Temp, shunt)

% Thermal voltage at room temperature in V
Vth = 0.02569;

% Temperature scaling factors
tmpJsc = (1 + tcJSC * (Temp - 25));
tmpVoc = (1 + tcVOC * (Temp - 25));

% JSC for cell at temperature of the cell Temp
jsc = jsc_RT .* tmpJsc;

% Voc
thrs = jsc > 100*j0;
if strcmp(shunt,'without') % w/o shunt
    VOC_RT = n * Vth * log(jsc./j0);
elseif strcmp(shunt,'with')  % with shunt
    VOC_RT = jsc/1000 * Rsh - n * Vth * lambertwlog( log(j0/1000 * Rsh) + Rsh * (jsc + j0)/(1000 * n * Vth) - log((n * Vth)) ) + j0/1000 * Rsh;
end
VOC = thrs .* VOC_RT .* tmpVoc;
VOC(isnan(VOC))=0;

% Define V
V = linspace(-0.2,VOC*1.1,100);

% Calculate J
if strcmp(shunt,'without') % w/o shunt
    z = (Rs*j0/1000)/(n*Vth) * exp( (Rs*(jsc+j0)/1000+V-VOC+VOC_RT) / (n*Vth) );
    j = 1000 * ( (jsc + j0)/1000 - n*Vth/Rs * lambertw(z) );
elseif strcmp(shunt,'with')  % with shunt
    z = (Rs*j0/1000)/(n*Vth*(1+Rs/Rsh)) * exp( (Rs*(jsc+j0)/1000+V-VOC+VOC_RT) / (n*Vth*(1+Rs/Rsh)) );
    j = 1000 * ( ( jsc/1000 + j0/1000 -(V-VOC+VOC_RT)/Rsh)/(1+Rs/Rsh) - n*Vth/Rs * lambertw(z) );
end

% Power
P = j .* V;
[P_max,I] = max(P,[],2);

% Generate output
JMPP = j(I)';%[j] = mA/cm^2
for z=1:length(I)
    VMPP(z,:) = V(z,I(z));%[V] = V
end
FF = JMPP .* VMPP ./ (jsc .* VOC);%FF dimensionless
FF(isnan(FF))=0;
P_el = P_max.*10;%[P] = W/m^2

end