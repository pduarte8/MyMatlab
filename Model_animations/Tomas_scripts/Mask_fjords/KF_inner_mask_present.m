% Create a mask layer for inner part of Kongsfjorden, excluding areas
% in front of Kronebreen and Kongsbreen tidewater glaciers, along
% transects across Loven Islands, Dyrevika, Kongsbreen North and 
% Kronebreen, based on existing land mask and a polygon covering the
% Kongsfjorden area. 
% The limits of the transects are defined by the lines:
%
%   South coast near Isgrunnen   (78.8906 N, 12.3452 E)
%   Blomstrand peninsula (south) (78.9725 N, 12.1876 E)
%   Blomstrand peninsula (north) (78.9953 N, 12.2001 E)
%   North coast (west)           (79.0025 N, 12.2495 E)
%   North coast (east)           (78.9932 N, 12.4709 E)
%   Ossian Sars mountain (north) (78.9704 N, 12.4657 E)
%   Ossian Sars mountain (south) (78.9139 N, 12.4847 E)
%   South coast (east)           (78.8780 N, 12.4056 E)

clear;

% Grid file
gname = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/kongsfjorden_160m_grid_present.nc';

% Get grid dimensions
xi = length(ncread(gname,'xi_rho'));
eta = length(ncread(gname,'eta_rho'));

% Get grid
lat = ncread(gname,'lat_rho');
lon = ncread(gname,'lon_rho');
mask = ncread(gname,'mask_rho');

% Inner polygon in front of Kongsvegen and Kronebreen

%lat_poly = [  78.892  78.945  78.929  78.844  78.892  ];
%lon_poly = [  12.180  12.460  12.805  12.469  12.180 ];

% Inner polygon in front of Kongsbreen and Conwaybreen

lat_poly = [  78.945  78.929  79.07  78.945 ];
lon_poly = [  12.460  12.805  12.366 12.460 ];


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
save('mask_KF_inner_present.mat', 'mask_fjord');
