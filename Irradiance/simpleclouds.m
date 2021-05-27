% ------------- DOKUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:          Derives the real spectra considering cloud
%                        coverage and real diffuse spectra
%
% #INPUT:                Code_location,Alias_location
%
% #OUTPUT:               /
%
% #SAVED DATA:           Saves clear Sky Irradiance 'Irr_spectra_clouds.mat'.
%                         FORMAT: Saved as *.mat file containing three matrixes (365x24, length([280:1:4000)).  
%                         1st matrix: 'Irr_spectra_clouds_wavelength'
%                         2nd matrix: 'Irr_spectra_clouds_direct_horizontal'
%                         3rd matrix: 'Irr_spectra_clouds_diffuse_horizontal'
%
% #REQUIRED SUBFUNCTIONS: /
%
% #ADD COMMENTS:          /
%
% -----------------------------------------------------------
%
function simpleclouds(Code_location,Alias_location)

NAME_location = [Code_location,'_',Alias_location];
folder_target_data = ['Irradiance\Spectra_',NAME_location];
load([folder_target_data,'\TMY3_',NAME_location,'.mat']);
load([folder_target_data,'\Irr_spectra_clear_sky.mat']);
clear Irr_spectra_clear_sky_direct_normal

Irr_spectra_clouds_wavelength         = Irr_spectra_clear_sky_wavelength;                                    
Irr_spectra_clouds_direct_horizontal  = zeros(size(Irr_spectra_clear_sky_diffuse_horizontal));  
Irr_spectra_clouds_diffuse_horizontal = zeros(size(Irr_spectra_clear_sky_diffuse_horizontal));


for index2=1:size(Irr_spectra_clear_sky_diffuse_horizontal,1)

  % Calculate impact of clouds 
  if sum(Irr_spectra_clear_sky_direct_horizontal(index2,:))> 0
      % Speculare Irradiance Spectra as "direct horizontal irradiance spectra"     
      Irr_spectra_clouds_direct_horizontal(index2,:) = ...
          Irr_spectra_clear_sky_direct_horizontal(index2,:)./sum(Irr_spectra_clear_sky_direct_horizontal(index2,:)).*...
          max(0,(Data_TMY3(index2,11)-Data_TMY3(index2,13)));       % max because could be negative
     
      if sum(Irr_spectra_clear_sky_diffuse_horizontal(index2,:))> 0
      % Diffuse Irradiance Spectra as "diffuse horizontal irradiance spectra"
      Irr_spectra_clouds_diffuse_horizontal(index2,:) = ...
          ((Irr_spectra_clear_sky_diffuse_horizontal(index2,:).*((10-Data_TMY3(index2,23))/10) + Irr_spectra_clear_sky_direct_horizontal(index2,:).*(Data_TMY3(index2,23))/10)) ... 
          ./sum(((Irr_spectra_clear_sky_diffuse_horizontal(index2,:).*((10-Data_TMY3(index2,23))/10) + Irr_spectra_clear_sky_direct_horizontal(index2,:).*((Data_TMY3(index2,23))/10))))...
          .*((Data_TMY3(index2,13)));
      end
      
  else
      % Speculare Irradiance Spectra as "direct horizontal irradiance spectra"     
      Irr_spectra_clear_sky_direct_horizontal(index2,:) = zeros(size(Irr_spectra_clear_sky_direct_horizontal(index2,:)));
  
      if sum(Irr_spectra_clear_sky_diffuse_horizontal(index2,:))> 0
      % Diffuse Irradiance Spectra as "diffuse horizontal irradiance spectra"
      Irr_spectra_clouds_diffuse_horizontal(index2,:) = ...
          ((Irr_spectra_clear_sky_diffuse_horizontal(index2,:).*((10-Data_TMY3(index2,23))/10) + Irr_spectra_clear_sky_direct_horizontal(index2,:).*(Data_TMY3(index2,23))/10)) ... 
          ./sum(((Irr_spectra_clear_sky_diffuse_horizontal(index2,:).*((10-Data_TMY3(index2,23))/10) + Irr_spectra_clear_sky_direct_horizontal(index2,:).*((Data_TMY3(index2,23))/10))))...
          .*((Data_TMY3(index2,13)));
      end
  end
  
%         figure
%         hold on
%         plot(Irr_spectra_clouds_wavelength(index2,:), Irr_spectra_clouds_diffuse_horizontal(index2,:))
%         plot(Irr_spectra_clouds_wavelength(index2,:), Irr_spectra_clouds_direct_horizontal(index2,:))
%         hold off
end



save([folder_target_data,'\Irr_spectra_clouds.mat'], 'Irr_spectra_clouds_wavelength','Irr_spectra_clouds_direct_horizontal','Irr_spectra_clouds_diffuse_horizontal');


end

