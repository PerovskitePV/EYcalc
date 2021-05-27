% ------------- DOKUMENTATION OF THIS FUNCTION -------------
%         
% #DESCRIPTION:           Calculates the Fresnel Reflection and Transmission at an interface between 
%                         a medium with refractive index n1 and a medium with refractive index n2 for
%                         all angles of incidence and the whole spectral range of n1 and n2. The 
%                         polarization is taken to be pol (can be 's' or 'p', or anything for 
%                         unpolarized light). n1 and n2 must be vectors.
%
% #INPUT:                 n1 - refractive index of first layer
%                         n2 - refractive index of second layer
%                         pol - polarization (TE, TM or mixed)
%                         lambda - wavelength vector
%
%
% #OUTPUT:                R,T with size [90 x lambda]
%
% #SAVED DATA:            -
%
% #REQUIRED SUBFUNCTIONS: -
%
% #ADD COMMENTS:          If we start in absorbing media: isreal(n1)=false, we use TMM to calculate proper R and T
%
% -----------------------------------------------------------

function [R,T] = Fresnel(n1,n2,pol,lambda)
    
    % If n1=n2, T=1 and R=0
    if sum(n1-n2)==0
        T = ones(90,length(n1));
        R = zeros(90,length(n1));

    % If n1 and n2 are different, calculate the Fresnel R,T matrix
    % First, check if n1 is real:
    else%if sum(n1-n2)~=0 && isreal(n1)
        phi1 = 0:89;
        phi2 = real(floor( asind( sind(phi1)'*(n1./n2) )));
        % handel TIR
        phi2(phi2>=90)=NaN;

        rs = (cosd(phi1)'*n1 - n2.*cosd(phi2))./(cosd(phi1)'*n1+n2.*cosd(phi2));
        rp = (cosd(phi1)'*n2 - n1.*cosd(phi2))./(cosd(phi1)'*n2+n1.*cosd(phi2));

        switch pol
            case 'TE'
                R = abs(rs).^2;
            case 'TM'
                R = abs(rp).^2;
            otherwise
                R = ((abs(rs).^2 + abs(rp).^2) * .5 );
        end

        R(isnan(R))=1;
        T = 1-R;

%     % If n1 is complex, we use the TMM code to calculate the R and T properties    
%     else
%         phi1 = 0:89;
%         [R_TMM_TE, T_TMM_TE, R_TMM_TM, T_TMM_TM] = deal(zeros(length(phi1),length(n1)));
%         
%         for aoi=1:length(phi1)   
%             [R_TMM_TE(aoi,:),~,T_TMM_TE(aoi,:),~] = coh_tmm('TE', [n1.' n2.' n2.'], [Inf 10 Inf], phi1(aoi),lambda);
%             [R_TMM_TM(aoi,:),~,T_TMM_TM(aoi,:),~] = coh_tmm('TM', [n1.' n2.' n2.'], [Inf 10 Inf], phi1(aoi),lambda);
%         end
%         
%         R = ( R_TMM_TE + R_TMM_TM ) /2;
%         T = ( T_TMM_TE + T_TMM_TM ) /2;
    end


end