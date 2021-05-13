% Compute surface area and volume for Kongsfjorden subregion

clear all;

inPath = check_path('D:\Data\ROMS\matlab\utility');
if (~inPath)
    disp('Add ROMS toolbox to path');
    run('D:\Data\ROMS\matlab\startup.m');
end

% Grid file
gname = 'D:\Data\ROMS\norkyst\NorFjords-160m_Forcing\Kongsfjorden-160m_v2_future\Grid\kongsfjorden_160m_grid_future.nc';

% Fjord subregion mask
load('mask_KF_Blomstrand_future.mat');

% Get grid dimensions
xi = length(ncread(gname,'xi_rho'));
eta = length(ncread(gname,'eta_rho'));

% Get grid
lat_u = ncread(gname,'lat_u');
lon_u = ncread(gname,'lon_u');
lat_v = ncread(gname,'lat_v');
lon_v = ncread(gname,'lon_v');

% Get depth
h = ncread(gname,'h');

KF_Area = 0;
KF_Volume = 0;
for i=1:xi
    for j=1:eta
        if mask_fjord(i,j)   % true if 1
            dx = geodesic_dist(lon_u(i,j),lat_u(i,j),lon_u(i,j+1),lat_u(i,j+1),1);
            dy = geodesic_dist(lon_v(i,j),lat_v(i,j),lon_v(i+1,j),lat_v(i+1,j),1);
            KF_Area = KF_Area + dx*dy;
            KF_Volume = KF_Volume + h(i,j)*dx*dy;
        end
    end
end

save('volar_KF_Blomstrand_future.mat', 'KF_Area', 'KF_Volume');

KF_Area = KF_Area.*1.E-6;
KF_Volume = KF_Volume.*1.E-9;

disp(sprintf('Area = %.2f km^2 ; Volume = %.2f km^3', KF_Area, KF_Volume));