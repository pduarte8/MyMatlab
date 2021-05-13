% Create a mask layer for outer part of Kongsfjorden, along transects
% from Kvadehuken to Kapp Guissez, across Blomstrandhamna, and from
% Blomstrand peninsula to Ny-Alesund, based on existing land mask and 
% a polygon covering the Kongsfjorden area. 
% The limits of the transects are defined by the lines:
%
%   Kvadehuken (south coast)     (78.9698 N, 11.3677 E)
%   Kapp Guissez (north coast)   (79.0673 N, 11.6528 E)
%   North coast (east)           (79.0042 N, 12.0184 E)
%   Blomstrand peninsula (north) (78.9935 N, 12.0164 E)
%   Blomstrand peninsula (south) (78.963  N, 12.0267 E)
%   NyAlesund                    (78.9288 N, 11.9094 E)

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
lat_poly = [  78.92    78.9698  79.0673  79.0673  79.0042  ...
              78.9935  78.963   78.9288  78.92    78.92    ];
lon_poly = [  11.6     11.3677  11.6528  12.0184  12.0184  ...
              12.0164  12.0267  11.9094  11.9094  11.6     ];

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
save('mask_KF_outer_present.mat', 'mask_fjord');
