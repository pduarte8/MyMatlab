% Create a mask layer for middle part of Kongsfjorden south of 
% Blomstrand peninsula, along transects from Ny-Alesund to Blomstrand 
% peninsula and across Loven Islands, based on existing land mask and 
% a polygon covering the Kongsfjorden area. 
% The limits of the transects are defined by the lines:
%
%   NyAlesund                   (78.9288 N, 11.9094 E)
%   Blomstrand peninsula (west) (78.963  N, 12.0267 E)
%   Blomstrand peninsula (east) (78.9725 N, 12.1876 E)
%   South coast near Isgrunnen  (78.8906 N, 12.3452 E)

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
lat_poly = [  78.88    78.9288  78.963   78.98    78.98  ...
              78.9725  78.8906  78.88    78.88    ];
lon_poly = [  11.9094  11.9094  12.0267  12.0267  12.1876  ...
              12.1876  12.3452  12.3452  11.9094  ];

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
save('mask_KF_middle_south_present.mat', 'mask_fjord');
