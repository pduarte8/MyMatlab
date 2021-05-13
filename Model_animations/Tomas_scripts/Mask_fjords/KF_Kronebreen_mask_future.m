% Create a mask layer for area in front of Kronebreen until 
% Kronebreen transect,
% based on existing land mask and a % polygon covering the Kongsfjorden 
% area. The limits of the transect is defined by the line:
%
%   South coast  (78.8780 N, 12.4056 E)
%   North coast  (78.9139 N, 12.4847 E)

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
lat_poly = [lat(1,1) lat(120,1) 78.90 78.90 78.93 78.95  ...
            78.9139 78.8780                              ...  % transect
            78.8 lat(1,1)];
lon_poly = [lon(1,1) lon(120,1) 12.90 12.67 12.64 12.50  ...
            12.4847 12.4056                              ...  % transect
            12.4 lon(1,1)];
lat_poly = [  78.879487  78.913930  78.929000  78.751000   78.879487];
lon_poly = [  12.405971  12.476839  13.053000  13.004000   12.405971];
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
save('mask_KF_Kronebreen_future.mat', 'mask_fjord');
