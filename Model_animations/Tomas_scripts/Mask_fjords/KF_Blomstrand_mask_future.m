% Create a mask layer for inner part of Kongsfjorden, along transects
% from Ny-Alesund to Blomstrand peninsula and across Blomstrand bay,
% based on existing land mask and a % polygon covering the Kongsfjorden 
% area. The limits of the transects are defined by the lines:
%
%   NyAlesund            (78.9288 N, 11.9094 E)
%   Blomstrand peninsula (78.963  N, 12.0267 E)
%   Blomstrand bay south (78.9935 N, 12.0164 E)
%   Blomstrand bay north (79.0042 N, 12.0184 E)

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
lat_poly = [lat(1,1) lat(212,1) 79.06        ...
            79.0042 78.9935 78.963  78.9288  ...
            78.88 lat(1,1)];
lon_poly = [lon(1,1) lon(212,1) 12.00        ...
            12.0184 12.0164 12.0267 11.9094  ...
            11.9 lon(1,1)];

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
save('mask_KF_Blomstrand_future.mat', 'mask_fjord');
