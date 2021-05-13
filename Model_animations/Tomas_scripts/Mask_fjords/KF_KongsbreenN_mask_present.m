% Create a mask layer for area in front of Kongsbreen North until 
% Kongsbreen transect,
% based on existing land mask and a % polygon covering the Kongsfjorden 
% area. The limits of the transect is defined by the line:
%
%   South coast  (78.9704 N, 12.4657 E)
%   North coast  (78.9932 N, 12.4709 E)

clear;

% Grid file
gname = 'D:\Data\ROMS\norkyst\NorFjords-160m_Forcing\Kongsfjorden-160m_v2_present\Grid\kongsfjorden_160m_grid_present_NEW.nc';

% Get grid dimensions
xi = length(ncread(gname,'xi_rho'));
eta = length(ncread(gname,'eta_rho'));

% Get grid
lat = ncread(gname,'lat_rho');
lon = ncread(gname,'lon_rho');
mask = ncread(gname,'mask_rho');

% Kongsfjorden polygon
lat_poly = [lat(120,1) lat(212,1) 79.05  ...
            78.9932 78.9704              ...  % transect
            78.95 78.93 78.90 78.90 lat(120,1)];
lon_poly = [lon(120,1) lon(212,1) 12.50  ...
            12.4709 12.4657              ...  % transect
            12.50 12.64 12.67 12.90 lon(120,1)];

% reshape lat & lon matrices
lat_vec = reshape(lat,[1,xi*eta]);
lon_vec = reshape(lon,[1,xi*eta]);

% create polygon mask
mask_vec = double(inpolygon(lat_vec,lon_vec,lat_poly,lon_poly));
mask_poly = reshape(mask_vec,[xi,eta]);

% create fjord mask
mask_fjord = mask .* mask_poly;

% plot mask
figure
subplot(2,1,1)
contourf(lon,lat,mask)
line(lon_poly,lat_poly,'Color','red')
%
subplot(2,1,2)
contourf(lon,lat,mask_fjord)

% save mask
save('mask_KF_KongsbreenN_present.mat', 'mask_fjord');
