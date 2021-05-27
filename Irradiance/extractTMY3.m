% ------------- DOKUMENTATION OF THIS FUNCTION -------------
% #DESCRIPTION:          Extracts data from TMY3 dataset and saves in Data_TMY3.mat
%
% #INPUT:                Code_location,Alias_location
%
% #OUTPUT:               /
%
% #SAVED DATA:           Relevant Data from the TMY3 data set.
%             FORMAT: Saved as *.mat file containing a Matrix called "Data_TMY3" in folder ['Spectra_'Code_location'_'Alias_Location'].  
%             ROWS:    365*24  - hourly resolved dataset - 365*24. 
%             COLUMNS: 22      - Defionition below
%             Data_TMY3(:,1)  - year   XXXX 
%             Data_TMY3(:,2)  - month  1-12
%             Data_TMY3(:,3)  - day    1-31
%             Data_TMY3(:,4)  - hour   0-23 
%             Data_TMY3(:,5)  - min    0-59   currently not relevant
%             Data_TMY3(:,6)  - sec    0-59   currently not relevant
%             Data_TMY3(:,7)  - angle to zenit, 0-180°, 0° if the sun is pointing towards the sky
%             Data_TMY3(:,8)  - azimut angle, 0-360°, 0° at north, 90° at east, 180° south ... 
%             Data_TMY3(:,9)  - Hourly extraterrestrial radiation normal to the sun  Wh/m2
%             Data_TMY3(:,10) - Hourly extraterrestrial radiation normal to the sun  Wh/m2
%             Data_TMY3(:,11) - Total global horizontal irradiance Wh/m2
%             Data_TMY3(:,12) - Total direct normal irradiance
%             Data_TMY3(:,13) - Total diffuse horizontal irradiance Wh/m2
%             Data_TMY3(:,14) - Dry-bulb temperature  C
%             Data_TMY3(:,15) - relative humuduty  %                                
%             Data_TMY3(:,16) - Station pressure mbar
%             Data_TMY3(:,17) - preciptable water 1/cm
%             Data_TMY3(:,18) - aoreosol optical depth [unitless]   used for turbidity 
%             Data_TMY3(:,19) - albedo [unitless] 
%             Data_TMY3(;,20)- longitude  
%             Data_TMY3(:,21)- latitude
%             Data_TMY3(:,22)- altitude 
%             Data_TMY3(:,23)- cloud coverage in tenth per cent   
%
% #REQUIRED SUBFUNCTIONS: 
%         1. sun_position(time, location); Code taken from NREL
%         2. lineArray = read_mixed_csv(fileName,delimiter)
%
% #ADD COMMENTS:    /
%
%-----------------------------------------------------------
%
function extractTMY3(Code_location,Alias_location)

NAME_location = [Code_location,'_',Alias_location];

input = read_mixed_csv(['Irradiance\Dataset_TMY3\',Code_location,'.CSV'],',/:');
location.longitude = str2double(input(1,6));
location.latitude  = str2double(input(1,5));
location.altitude  = str2double(input(1,7));

Data_TMY3 = zeros(size(input,1)-2,23);


% Note the colum number is +3 as we split the first two columns of the CSV
% file in 5 columns
%------------- Go 
Data_TMY3(:,9)  =  str2double(input(3:size(input,1),6));  % Hourly extraterrestrial radiation normal to the sun  Wh/m2
Data_TMY3(:,10) =  str2double(input(3:size(input,1),7));  % Hourly extraterrestrial radiation normal to the sun  Wh/m2
Data_TMY3(:,11) =  str2double(input(3:size(input,1),8));  % Global horizontal irradiance Wh/m2
Data_TMY3(:,12) =  str2double(input(3:size(input,1),11)); % Direct normal irradiance
Data_TMY3(:,13) =  str2double(input(3:size(input,1),14)); % Diffuse horizontal irradiance Wh/m2
Data_TMY3(:,14) =  str2double(input(3:size(input,1),35)); % Dry-bulb temperature  C
Data_TMY3(:,15) =  str2double(input(3:size(input,1),41)); % relative humuduty  %                                
Data_TMY3(:,16) =  str2double(input(3:size(input,1),44)); % Station pressure mbar
Data_TMY3(:,17) =  str2double(input(3:size(input,1),59)); % preciptable water 1/cm
Data_TMY3(:,18) =  str2double(input(3:size(input,1),62)); % aoreosol optical depth [unitless]   used for turbidity 
Data_TMY3(:,19) =  str2double(input(3:size(input,1),65)); % albedo [unitless]   
Data_TMY3(:,23) =  str2double(input(3:size(input,1),32)); %  cloud coverage in tenth per cent   
%------------- Time frame
 Data_TMY3(:,1) =  str2double(input(3:size(input,1),3));  % year   XXXX 
 Data_TMY3(:,2) =  str2double(input(3:size(input,1),1));  % month  1-12
 Data_TMY3(:,3) =  str2double(input(3:size(input,1),2));  % day    1-31 
 Data_TMY3(:,4) =  str2double(input(3:size(input,1),4));  % hour   0-23 
 Data_TMY3(:,5) =  str2double(input(3:size(input,1),5));  % min    0-59   currently not relevant
 Data_TMY3(:,6) =  str2double(input(3:size(input,1),5));  % sec    0-59   currently not relevant Just zeros !!!
 
for ii=1:1:(size(Data_TMY3,1))
 
 time.UTC = str2double(input(1,4));  
 
 time.year  =  Data_TMY3(ii,1);
 time.month =  Data_TMY3(ii,2);
 time.day   =  Data_TMY3(ii,3);
 time.hour  =  Data_TMY3(ii,4);
 time.min   =  Data_TMY3(ii,5);
 time.sec   =  Data_TMY3(ii,6);
 
 sun = sun_position(time, location);
 
 Data_TMY3(ii,7)  =  sun.zenith;          % angle to zenit, 0-180°, 0° if the sun is pointing towards the sky
 Data_TMY3(ii,8)  =  sun.azimuth;         % azimut angle, 0-360°, 0° at north, 90° at east, 180° south ... 
 Data_TMY3(ii,20) =  location.longitude;  % longitude  
 Data_TMY3(ii,21) =  location.latitude;   % latitude
 Data_TMY3(ii,22) =  location.altitude;   % altitude 
 
end

folder_target_data = ['Irradiance\Spectra_',NAME_location];
mkdir(folder_target_data);

% Save the data
save([folder_target_data,'\TMY3_', NAME_location,'.mat'], 'Data_TMY3');


end

