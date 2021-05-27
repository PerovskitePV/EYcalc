% ------------- DOKUMENTATION OF THIS FUNCTION -------------
% %DESCRIPTION:           Executes the calculation of the Poynting vector and the
%                         overall scattering matrix for all wavelengths and calculates reflection,
%                         absorption and transmission.
% %INPUT:                 n_array: CRI of the layers (array)
%                         d_list: Thicknesses of the layers (vector)
%                         pol: Polarisation of the field to be computed (string: 'TE' or 'TM')
%                         lam_vac: Wavelength for which the computation shall be conducted (vector)
%                         th_0: Angle of incidendence in degree
% %OUTPUT:                data: detailed overview of the simulation (struct)                       
% %SAVED DATA:            -
% %REQUIRED SUBFUNCTIONS: 
%
% %ADD COMMENTS:          -
% -----------------------------------------------------------
% Code is a modified version of:  https://github.com/hugadams/PAME/blob/master/pame/tmm_mod.py
%                                 https://github.com/sbyrnes321/tmm/blob/master/tmm_core.py
% -- Original License --
% Copyright (C) 2012-2020 Steven Byrnes
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial
% portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
% LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% -- 
%
% All credits go to Steven Byrnes!
% Adaptional and some minor changes have been done by Malte Langenhorst

function [R, A, T, th_f] = coh_tmm(pol, n_array, d_list, th_0, lam_vac)
    % """
    % Main "coherent transfer matrix method" calc. Given parameters of a stack,
    % calculates everything you could ever want to know about how light
    % propagates in it. (If performance is an issue, you can delete some of the
    % calculations without affecting the rest.)
    % 
    % pol is light polarization, "TE" or "TM".
    % 
    % n_list is the list of refractive indices, in the order that the light would
    % pass through them. The 0'th element of the list should be the semi-infinite
    % medium from which the light enters, the last element should be the semi-
    % infinite medium to which the light exits (if any exits).
    % 
    % th_0 is the angle of incidence: 0 for normal, pi/2 for glancing.
    % Remember, for a dissipative incoming medium (n_list[0] is not real), th_0
    % should be complex so that n0 sind(th0) is real (intensity is constant as
    % a function of lateral position).
    % 
    % d_list is the list of layer thicknesses (front to back). Should correspond
    % one-to-one with elements of n_list. First and last elements should be "inf".
    % 
    % lam_vac is vacuum wavelength of the light.
    % 
    % Outputs the following as a dictionary (see manual for details)
    % 
    % * r--reflection amplitude
    % * t--transmission amplitude
    % * R--reflected wave power (as fraction of incident)
    % * T--transmitted wave power (as fraction of incident)
    % * power_entering--Power entering the first layer, usually (but not always)
    % equal to 1-R (see manual).
    % * vw_list-- n'th element is [v_n,w_n], the forward- and backward-traveling
    % amplitudes, respectively, in the n'th medium just after interface with
    % (n-1)st medium.
    % * kz_list--normal component of complex angular wavenumber for
    % forward-traveling wave in each layer.
    % * th_list--(complex) propagation angle (in radians) in each layer
    % * pol, n_list, d_list, th_0, lam_vac--same as input
    % 
    % """
    % Convert lists to numpy arrays if they're not already.
    
    % Input tests
    if (size(th_0,2) > 1)
        error('This function is not vectorized; you need to run one calculation for one angle at a time');
    end
    if (size(n_array,1) ~= length(lam_vac)) || (size(d_list,1) ~= 1) || (size(n_array,2) ~= size(d_list,2))
        error("Problem with n_list or d_list!");
    end
    if (d_list(1) ~= Inf || d_list(end) ~= Inf)
        error('d_list must start and end with inf!');
    end
    if min(is_forward_angle(n_array(:,1), th_0)) == 0
        error('Error in n0 or th0!');
    end
    num_layers = size(n_array,2);
    num_lambdas = size(lam_vac,2);
    n_array(:,1) = real(n_array(:,1));
        
    % th_list is a list with, for each layer, the angle that the light travels
    % through the layer. Computed with Snell's law. Note that the "angles" may be
    % complex!
    th_array = list_snell(n_array, th_0);

    % kz is the z-component of (complex) angular wavevector for forward-moving
    % wave. Positive imaginary part means decaying.
    kz_array = 2 * pi * n_array .* cosd(th_array) ./ lam_vac';

    % delta is the total phase accrued by traveling through a given layer.
    % Ignore warning about inf multiplication
    delta = kz_array .* d_list;

    % For a very opaque layer, reset delta to avoid divide-by-0 and similar
    % errors. The criterion imag(delta) > 35 corresponds to sindgle-pass
    % transmission < 1e-30 --- small enough that the exact value doesn't
    % matter.
    for i=2:num_layers-1
        if max(imag(delta(:,i))) > 35
            delta(imag(delta(:,i)) > 35, i) = real(delta(imag(delta(:,i)) > 35, i)) + 35i;
            opacity_warning = true;
%             if exist('opacity_warning','var')
%                 disp(strcat("Warning: Layers that are almost perfectly opaque ",...
%                       "are modified to be slightly transmissive, ",...
%                       "allowing 1 photon in 10^30 to pass through. This is ",...
%                       "for numerical stability."));
%             end
        end
    end

    % t_list[i,j] and r_list[i,j] are transmission and reflection amplitudes,
    % respectively, coming from i, going to j. Only need to calculate this when
    % j=i+1. (2D array is overkill but helps avoid confusion.)
    t_list = zeros(num_layers, num_layers, num_lambdas);
    r_list = zeros(num_layers, num_layers, num_lambdas);
    for i = 1:num_layers-1
        t_list(i,i+1,:) = interface_t(pol, n_array(:,i), n_array(:,i+1), th_array(:,i), th_array(:,i+1));
        r_list(i,i+1,:) = interface_r(pol, n_array(:,i), n_array(:,i+1), th_array(:,i), th_array(:,i+1));
    end
    % At the interface between the (n-1)st and nth material, let v_n be the
    % amplitude of the wave on the nth side heading forwards (away from the
    % boundary), and let w_n be the amplitude on the nth side heading backwards
    % (towards the boundary). Then (v_n,w_n) = M_n (v_{n+1},w_{n+1}). M_n is
    % M_list[n]. M_0 and M_{num_layers-1} are not defined.
    % My M is a bit different than Sernelius's, but Mtilde is the same.
    M_list = zeros(num_layers,2,2,num_lambdas);
    for i=2:num_layers-1
        M_list(i,1,1,:) = reshape((1/t_list(i,i+1,:)),num_lambdas,1) .* exp(-1j*delta(:,i));
        M_list(i,1,2,:) = reshape((1/t_list(i,i+1,:)) .* r_list(i,i+1,:),num_lambdas,1) .* exp(-1j*delta(:,i));
        M_list(i,2,1,:) = reshape((1/t_list(i,i+1,:)) .* r_list(i,i+1,:),num_lambdas,1) .* exp(1j*delta(:,i));
        M_list(i,2,2,:) = reshape((1/t_list(i,i+1,:)),num_lambdas,1) .* exp(1j*delta(:,i));
    end
    Mtilde = repmat(eye(2),1,1,num_lambdas);
    for i=2:num_layers-1
        Mtilde_temp = Mtilde;
        Mtilde(1,1,:) = Mtilde_temp(1,1,:) .* reshape(M_list(i,1,1,:),1,1,num_lambdas) + Mtilde_temp(1,2,:) .* reshape(M_list(i,2,1,:),1,1,num_lambdas);
        Mtilde(1,2,:) = Mtilde_temp(1,1,:) .* reshape(M_list(i,1,2,:),1,1,num_lambdas) + Mtilde_temp(1,2,:) .* reshape(M_list(i,2,2,:),1,1,num_lambdas);
        Mtilde(2,1,:) = Mtilde_temp(2,1,:) .* reshape(M_list(i,1,1,:),1,1,num_lambdas) + Mtilde_temp(2,2,:) .* reshape(M_list(i,2,1,:),1,1,num_lambdas);
        Mtilde(2,2,:) = Mtilde_temp(2,1,:) .* reshape(M_list(i,1,2,:),1,1,num_lambdas) + Mtilde_temp(2,2,:) .* reshape(M_list(i,2,2,:),1,1,num_lambdas);
    end
    Mtilde_temp = Mtilde;
    Mtilde(1,1,:) = reshape(1./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(1,1,:) + reshape(r_list(1,2,:)./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(2,1,:);
    Mtilde(1,2,:) = reshape(1./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(1,2,:) + reshape(r_list(1,2,:)./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(2,2,:);
    Mtilde(2,1,:) = reshape(r_list(1,2,:)./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(1,1,:) + reshape(1./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(2,1,:);
    Mtilde(2,2,:) = reshape(r_list(1,2,:)./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(1,2,:) + reshape(1./t_list(1,2,:),1,1,num_lambdas) .* Mtilde_temp(2,2,:);
    
    % Net complex transmission and reflection amplitudes
    r = reshape(Mtilde(2,1,:)./Mtilde(1,1,:),num_lambdas,1);
    t = reshape(1./Mtilde(1,1,:),num_lambdas,1);

    % vw_list[n] = [v_n, w_n]. v_0 and w_0 are undefined because the 0th medium
    % has no left interface.
    vw_array = zeros(num_layers, 2, num_lambdas);
    vw = [t,zeros(num_lambdas,1)];
    vw_array(end,:,:) = vw.';
    for i = fliplr(1:num_layers-1)
        vw_temp = vw;
        vw(:,1) = reshape(M_list(i,1,1,:),num_lambdas,1) .* vw_temp(:,1) + reshape(M_list(i,1,2,:),num_lambdas,1) .* vw_temp(:,2);
        vw(:,2) = reshape(M_list(i,2,1,:),num_lambdas,1) .* vw_temp(:,1) + reshape(M_list(i,2,2,:),num_lambdas,1) .* vw_temp(:,2);
        vw_array(i,:,:) = vw.';
    end
    % Net transmitted and reflected power, as a proportion of the incoming light
    % power.
    R = R_from_r(r);
    T = T_from_t(pol, t, n_array(:,1), n_array(:,end), th_0, th_array(:,end));
    power_entering = power_entering_from_r(pol, r, n_array(:,1), th_0);
    th_f = round(th_array(:,end));
    if th_f == 91
        th_f = 90;
    end
    coh_tmm_data.r = r;
    coh_tmm_data.t = t;
    coh_tmm_data.R = R;
    coh_tmm_data.T = T;
    coh_tmm_data.power_entering = power_entering;
    coh_tmm_data.vw_array = vw_array;
    coh_tmm_data.kz_array = kz_array;
    coh_tmm_data.th_array = th_array;
    coh_tmm_data.pol = pol;
    coh_tmm_data.n_array = n_array;
    coh_tmm_data.d_list = d_list;
    coh_tmm_data.th_0 = th_0;
    coh_tmm_data.lam_vac = lam_vac;

%     dx = 0.01; % resolution of the calculation of poynting vector
%     ds = -10:dx:sum(d_list(2:end-1)); % #position in structure
%     coh_tmm_data.poyn = zeros(length(lam_vac),length(ds));
%     coh_tmm_data.absor = zeros(length(lam_vac),length(ds));
%     for d = 1:length(ds)
%         [layer, d_in_layer] = find_in_structure_with_inf(d_list, ds(d));
%         abs_tmm_data = position_resolved(layer, d_in_layer, coh_tmm_data);
%         coh_tmm_data.poyn(:,d) = coh_tmm_data.poyn;
%         coh_tmm_data.absor(:,d) = coh_tmm_data.absor;
%     end
    RAT = absorp_in_each_layer(coh_tmm_data); % first entry is reflection, last entry is transmission
    A = RAT(:,2:end-1)';
end


function angles = list_snell(n_list, th_0)
    % """
    % return list of angle theta in each layer based on angle th_0 in layer 0,
    % usindg Snell's law. n_list is index of refraction of each layer. Note that
    % "angles" may be complex!!
    % """
    % Important that the arcsind here is scipy.arcsind, not numpy.arcsind! (They
    % give different results e.g. for arcsind(2).)
    angles = asind(n_list(:,1)*sind(th_0) ./ n_list);
    % The first and last entry need to be the forward angle (the intermediate
    % layers don't matter, see https://arxiv.org/abs/1603.02720 Section 5)
    idx_first = is_forward_angle(n_list(:,1), angles(:,1));
    idx_last = is_forward_angle(n_list(:,end), angles(:,end));
    angles(idx_first==0,1) = 90 - angles(idx_first==0,1);
    angles(idx_last==0,end) = 90 - angles(idx_last==0,end);
end


function r = interface_r(pol, n_i, n_f, th_i, th_f)
    % """
    % reflection amplitude (from Fresnel equations)
    % 
    % polarization is either "s" or "p" for polarization
    % 
    % n_i, n_f are (complex) refractive index for incident and final
    % 
    % th_i, th_f are (complex) propegation angle for incident and final
    % (in radians, where 0=normal). "th" stands for "theta".
    % """
    if strcmp(pol,'TE')
        r = ((n_i .* cosd(th_i) - n_f .* cosd(th_f)) ./ (n_i .* cosd(th_i) + n_f .* cosd(th_f)));
    elseif strcmp(pol,'TM')
        r = ((n_f .* cosd(th_i) - n_i .* cosd(th_f)) ./ (n_f .* cosd(th_i) + n_i .* cosd(th_f)));
    else
        error("Polarization must be 'TE' or 'TM'")
    end
end


function t = interface_t(pol, n_i, n_f, th_i, th_f)
    % """
    % transmission amplitude (frem Fresnel equations)
    % 
    % polarization is either "s" or "p" for polarization
    % 
    % n_i, n_f are (complex) refractive index for incident and final
    % 
    % th_i, th_f are (complex) propegation angle for incident and final
    % (in radians, where 0=normal). "th" stands for "theta".
    % """
    if strcmp(pol,'TE')
        t = 2 * n_i .* cosd(th_i) ./ (n_i .* cosd(th_i) + n_f .* cosd(th_f));
    elseif strcmp(pol,'TM')
        t = 2 * n_i .* cosd(th_i) ./ (n_f .* cosd(th_i) + n_i .* cosd(th_f));
    else
        error("Polarization must be 'TE' or 'TM'");
    end
end


function R = R_from_r(r)
    % """
    % Calculate reflected power R, starting with reflection amplitude r.
    % """
    R = abs(r.*r);
end


function T = T_from_t(pol, t, n_i, n_f, th_i, th_f)
    % """
    % Calculate transmitted power T, starting with transmission amplitude t.
    % 
    % n_i,n_f are refractive indices of incident and final medium.
    % 
    % th_i, th_f are (complex) propegation angles through incident & final medium
    % (in radians, where 0=normal). "th" stands for "theta".
    % 
    % In the case that n_i,n_f,th_i,th_f are real, formulas simplify to
    % T=|t|^2 * (n_f cosd(th_f)) / (n_i cosd(th_i)).
    % 
    % See manual for discussion of formulas
    % """
    if strcmp(pol,'TE')
        T = abs(t.*t) .* ((real(n_f.*cosd(th_f))) ./ real(n_i*cosd(th_i)));
    elseif strcmp(pol,'TM')
        T = abs(t.*t) .* ((real(n_f.*conj(cosd(th_f)))) ./ real(n_i*conj(cosd(th_i))));
    else
        error("Polarization must be 'TE' or 'TM'");
    end
end


function final_answer = absorp_in_each_layer(coh_tmm_data)
    % """
    % An array listing what proportion of light is absorbed in each layer.
    % 
    % Assumes the final layer eventually absorbs all transmitted light.
    % 
    % Assumes the initial layer eventually absorbs all reflected light.
    % 
    % Entries of array should sum to 1.
    % 
    % coh_tmm_data is output of coh_tmm()
    % """
    num_lambdas = length(coh_tmm_data.lam_vac);
    num_layers = length(coh_tmm_data.d_list);
    power_entering_each_layer = zeros(num_lambdas,num_layers);
    power_entering_each_layer(:,1) = ones(num_lambdas,1);
    power_entering_each_layer(:,2) = coh_tmm_data.power_entering;
    power_entering_each_layer(:,end) = coh_tmm_data.T;
    for i=3:num_layers-1
        data = position_resolved(i, 0, coh_tmm_data);
        power_entering_each_layer(:,i) = data.poyn;
    end
    final_answer = zeros(num_lambdas,num_layers);
    final_answer(:,1:end-1) = -diff(power_entering_each_layer,1,2);
    final_answer(:,end) = power_entering_each_layer(:,end);
    final_answer(abs(final_answer) < eps) = 0;
end


function answer = is_forward_angle(n, theta)
    %     """
    %     if a wave is traveling at angle theta from normal in a medium with index n,
    %     calculate whether or not this is the forward-traveling wave (i.e., the one
    %     going from front to back of the stack, like the incoming or outgoing waves,
    %     but unlike the reflected wave). For real n & theta, the criterion is simply
    %     -pi/2 < theta < pi/2, but for complex n & theta, it's more complicated.
    %     See https://arxiv.org/abs/1603.02720 appendix D. If theta is the forward
    %     angle, then (pi-theta) is the backward angle and vice-versa.
    %     """
    if max(real(n) .* imag(n) < 0)
        error(strcat("For materials with gain, it's ambiguous which ", ...
                                  "beam is incoming vs outgoing. See ",...
                                  "https://arxiv.org/abs/1603.02720 Appendix C.\n", ...
                                  "n: ", num2str(n), "   angle: ", num2str(theta)));
    end
    if length(theta) == 1
        theta = ones(length(n),1) * theta;
    end
    ncostheta = n .* cosd(theta);
    answer = NaN(size(n));
    % Either evanescent decay or lossy medium. Either way, the one that
    % decays is the forward-moving wave
    answer(abs(imag(ncostheta)) > 100 * eps) = (imag(ncostheta(abs(imag(ncostheta)) > 100 * eps)) > 0);
    % Forward is the one with positive Poynting vector
    % Poynting vector is Re[n cosd(theta)] for s-polarization or
    % Re[n cosd(theta*)] for p-polarization, but it turns out they're consistent
    % so I'll just assume s then check both below
    answer(abs(imag(ncostheta)) <= 100 * eps) = (real(ncostheta(abs(imag(ncostheta)) <= 100 * eps)) > 0);
    % convert from numpy boolean to the normal Python boolean
    % double-check the answer ... can't be too careful!
    error_string = strcat("It's not clear which beam is incoming vs outgoing. Weird", ...
                    " index maybe?\n", ...
                    "n: ", num2str(n), "   angle: ", num2str(theta));
    if min(answer) == true
        assert(min(imag(ncostheta) > -100 * eps), error_string);
        assert(min(real(ncostheta) > -100 * eps), error_string);
        assert(min((n .* real(cosd(conj(theta)))) > -100 * eps), error_string);
    else
        assert(min(imag(ncostheta) < 100 * eps), error_string);
        assert(min(real(ncostheta) < 100 * eps), error_string);
        assert(min((n .* real(cosd(conj(theta)))) < 100 * eps), error_string);
    end
end


function data = position_resolved(layer, distance, coh_tmm_data)
    % """
    % Starting with output of coh_tmm(), calculate the Poynting vector,
    % absorbed energy density, and E-field at a specific location. The
    % location is defined by (layer, distance), defined the same way as in
    % find_in_structure_with_inf(...).
    % 
    % Returns a dictionary containing:
    % 
    % * poyn - the component of Poynting vector normal to the interfaces
    % * absor - the absorbed energy density at that point
    % * Ex and Ey and Ez - the electric field amplitudes, where
    % z is normal to the interfaces and the light rays are in the x,z plane.
    % 
    % The E-field is in units where the incoming |E|=1; see
    % https://arxiv.org/pdf/1603.02720.pdf for formulas.
    % """
    num_lambdas = length(coh_tmm_data.lam_vac);
    if layer > 1
        v = reshape(coh_tmm_data.vw_array(layer,1,:),num_lambdas,1);
        w = reshape(coh_tmm_data.vw_array(layer,2,:),num_lambdas,1);
    else
        v = ones(length(coh_tmm_data.lam_vac),1);
        w = coh_tmm_data.r;
    end
    kz = coh_tmm_data.kz_array(:,layer);
    th = coh_tmm_data.th_array(:,layer);
    n = coh_tmm_data.n_array(:,layer);
    n_0 = coh_tmm_data.n_array(:,1);
    th_0 = coh_tmm_data.th_0;
    pol = coh_tmm_data.pol;

    assert((layer >= 1 && 0 <= distance <= coh_tmm_data.d_list(layer)) || (layer == 0 && distance <= 0));

    % Amplitude of forward-moving wave is Ef, backwards is Eb
    Ef = v .* exp(1j * kz * distance);
    Eb = w .* exp(-1j * kz * distance);

    % Poynting vector
    if strcmp(pol,'TE')
        poyn = (real(n.*cosd(th).*conj(Ef+Eb).*(Ef-Eb))) ./ real(n_0.*cosd(th_0));
    elseif strcmp(pol,'TM')
        poyn = ((real(n.*conj(cosd(th)).*(Ef+Eb).*conj(Ef-Eb))) ./ real(n_0*conj(cosd(th_0))));
    end

    % Absorbed energy density
    if strcmp(pol,'TE')
        absor = imag(n.*cosd(th).*kz).*abs(Ef+Eb).^2 ./ real(n_0*cosd(th_0));
    elseif strcmp(pol,'TM')
        absor = imag(n.*conj(cosd(th)).*(kz.*abs(Ef-Eb).^2-conj(kz).*abs(Ef+Eb).^2)) ./ real(n_0*conj(cosd(th_0)));
    end
                
    % Electric field
    if strcmp(pol,'TE')
        Ex = zeros(num_lambdas,1);
        Ey = Ef + Eb;
        Ez = zeros(num_lambdas,1);
    elseif strcmp(pol,'TM')
        Ex = (Ef - Eb) .* cosd(th);
        Ey = zeros(num_lambdas,1);
        Ez = (-Ef - Eb) .* sind(th);
    end
    data.poyn = poyn;
    data.absor = absor;
    data.Ex = Ex;
    data.Ey = Ey;
    data.Ez = Ez;
end


function P_in = power_entering_from_r(pol, r, n_i, th_i)
    %     """
    %     Calculate the power entering the first interface of the stack, starting with
    %     reflection amplitude r. Normally this equals 1-R, but in the unusual case
    %     that n_i is not real, it can be a bit different than 1-R. See manual.
    % 
    %     n_i is refractive index of incident medium.
    % 
    %     th_i is (complex) propegation angle through incident medium
    %     (in radians, where 0=normal). "th" stands for "theta".
    %     """
    if strcmp(pol,'TE')
        P_in = (real(n_i*cosd(th_i).*(1+conj(r)).*(1-r)) ./ real(n_i*cosd(th_i)));
    elseif strcmp(pol,'TM')
        P_in = (real(n_i*conj(cosd(th_i)).*(1+r).*(1-conj(r))) ./ real(n_i*conj(cosd(th_i))));
    else
        error("Polarization must be 'TE' or 'TM'");
    end
end
    