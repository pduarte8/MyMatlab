% Create a mask layer for Kongsfjorden, based on existing land mask and a
% polygon covering the Kongsfjorden area. The limit of Kongsfjorden is
% defined by a line:
%   south coast (78.97 N, 11.37 E)
%   north coast (79.07 N, 11.66 E)

clear;

% Grid file
gname = 'D:\Data\ROMS\norkyst\NorFjords-160m_Forcing\Kongsfjorden-160m_v2_future\Grid\kongsfjorden_160m_grid_future.nc';

% Get grid dimensions
xi = length(ncread(gname,'xi_rho'));
eta = length(ncread(gname,'eta_rho'));

% Get grid
lat = ncread(gname,'lat_rho');
lon = ncread(gname,'lon_rho');
mask = ncread(gname,'mask_rho');

% Kongsfjorden polygon
lat_poly = [lat(1,1) lat(212,1) 79.07 78.97 lat(1,1)];
lon_poly = [lon(1,1) lon(212,1) 11.66 11.37 lon(1,1)];

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
save('mask_kongsfjorden_future.mat', 'mask_fjord');
