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
% #ADD COMMENTS: -    
% -----------------------------------------------------------

function [VOC,FF,P_el,JMPP,VMPP,j,V] = calcsingleelectrics(jsc_RT, Rs, Rsh, j0, n, tcJSC, tcVOC, Temp, shunt)

% Thermal voltage at room temperature in V
Vth = 0.02569;

% Temperature scaling factors
tmpJsc = (1 + tcJSC * (Temp - 25));
tmpVoc = (1 + tcVOC * (Temp - 25));

% JSC for cell at temperature of the cell Temp
jsc = jsc_RT .* tmpJsc;

% VOC for cell at room temperature & at cell temperatures
thrs = jsc > 100*j0;

if strcmp(shunt,'without')  % w/o shunt
    VOC_RT = n * Vth * log(jsc./j0);
elseif strcmp(shunt,'with')  % with shunt
    VOC_RT = jsc/1000 * Rsh - n * Vth * lambertwlog( log(j0/1000 * Rsh) + Rsh * (jsc + j0)/(1000 * n * Vth) - log((n * Vth)) ) + j0/1000 * Rsh;
end
VOC = thrs .* VOC_RT .* tmpVoc;
VOC(isnan(VOC))=0;

% Current density in mA/cm^2
j = linspace(-0.2,1.05*max(jsc),100);

% Calculate the voltage points
if strcmp(shunt,'without')  % w/o shunt
    Vcond = ( repmat(jsc, [1 length(j)]) - repmat(j,[length(jsc) 1]) )/j0;
    V = n * Vth * log((Vcond>0).*Vcond) - j/1000 * Rs - VOC_RT + VOC;
    V(Vcond<=0)=0;
elseif strcmp(shunt,'with')  % with shunt
    lambw = log(j0 * Rsh) - log(1000 * n * Vth) + Rsh * (-j + jsc + j0)/(1000 * n * Vth);
    V = - j/1000 * ( Rs + Rsh) + jsc/1000 * Rsh - n * Vth .* lambertwlog( lambw ) +  j0/1000 * Rsh - VOC_RT + VOC;    
end

P = j .* V;
[P_max,I] = max(P,[],2); P_max(P_max<0)=0;

% Generate output
JMPP = j(I)'; JMPP(P_max==0) = 0;%[j] = mA/cm^2
for z=1:length(I)
    VMPP(z,:) = V(z,I(z));%[V] = V
end
VMPP(VMPP<0)=0;
FF = JMPP .* VMPP ./ (jsc .* VOC);%FF dimensionless
FF(isnan(FF))=0; FF(isinf(FF))=0; FF(FF<0);
P_el = P_max.*10;%[P] = W/m^2

end