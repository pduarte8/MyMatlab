% Create a mask layer for middle part of Kongsfjorden north of 
% Blomstrand peninsula, along transects across Blomstrandhamna
% and across Dyrevika, based on existing land mask and a
% polygon covering the Kongsfjorden area. 
% The limits of the transects are defined by the lines:
%
%   Blomstrand peninsula (west) (78.9935 N, 12.0164 E)
%   North coast (west)          (79.0042 N, 12.0184 E)
%   North coast (east)          (79.0025 N, 12.2495 E)
%   Blomstrand peninsula (east) (78.9953 N, 12.2001 E)

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
lat_poly = [  78.98    78.9935  79.0042  79.04    79.04    ...
              79.0025  78.9953  78.98    78.98    ];
lon_poly = [  12.0267  12.0164  12.0184  12.0184  12.2495  ...
              12.2495  12.2001  12.1876  12.0267  ];

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
save('mask_KF_middle_north_present.mat', 'mask_fjord');
