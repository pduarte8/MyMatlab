% Extract data for comparison with pyplume model
%
% Kongsvegen : index 347
% Kronebreen : index 351
% Kongsbreen : index 362
%
clear;

inPath = check_path('D:\Data\ROMS\matlab\utility');
if (~inPath)
    disp('Add ROMS toolbox to path');
    run('D:\Data\ROMS\matlab\startup.m');
end

A = load('D:\Data\ROMS\setup_Kongsfjord-K160\K160m_v2_present\Rivers\Liston_SnowModel_v2\Edge_Q_fluxes.dat');
B = load('D:\Data\ROMS\setup_Kongsfjord-K160\K160m_v2_present\Rivers\Liston_SnowModel_v2\Edge_Q_points_kongsfjorden_160m_present.txt');
%gname = 'D:\Data\ROMS\norkyst\NorFjords-160m_Forcing\Kongsfjorden-160m_v2_present\Grid\kongsfjorden_160m_grid_present_NEW.nc';
%gname = 'D:\Data\ROMS\setup_Kongsfjord-K160\K160m_v2_present\Grid\kongsfjorden_160m_grid_present.nc';
gname = 'D:\Data\ROMS\setup_Kongsfjord-K160\K160m_v3_present\Grid\kongsfjorden_160m_grid_present.nc';

rdir = 'D:\Data\roms_run\Kongsfjorden-160m\K160_v2\present_subglacial\daily_avg\'
prename = 'kongsfjorden_160m_avg.nc4_';
%nfiles = 30;
nfiles = 122;
t0 = datenum([2007 05 31 0 0 0]);  % Last day before beginning of data record

TEMP_1 = zeros(nfiles,35);
SALT_1 = zeros(nfiles,35);
TEMP_2 = zeros(nfiles,35);
SALT_2 = zeros(nfiles,35);
%TEMP_3 = zeros(nfiles,35);
%SALT_3 = zeros(nfiles,35);

% Extract dates and fluxes
time = datenum(A(:,1:3));
Q_flux = A(:,4:end);

% Extract Q for main glaciers
%Q_kongsvegen = Q_flux(:,347);
%Q_kronebreen = Q_flux(:,351);
%Q_kongsbreen = Q_flux(:,362);

i1 = find(time >= datenum([2007,6,1]),1,'first');
i2 = find(time < datenum([2007,10,1]),1,'last');
[y,m,d] = datevec(time(i1:i2));

% Extract Q for main glaciers
Q_kongsvegen = [y m d Q_flux(i1:i2,347)];
Q_kronebreen = [y m d Q_flux(i1:i2,351)];
Q_kongsbreen = [y m d Q_flux(i1:i2,362)];

% write discharge to file
fid = fopen('Kongsvegen_discharge.csv','w');
fprintf(fid,'year, month, day, m3/s\n');
fclose(fid);
dlmwrite('Kongsvegen_discharge.csv',Q_kongsvegen,'-append');

fid = fopen('Kronebreen_discharge.csv','w');
fprintf(fid,'year, month, day, m3/s\n');
fclose(fid);
dlmwrite('Kronebreen_discharge.csv',Q_kronebreen,'-append');

fid = fopen('Kongsbreen_discharge.csv','w');
fprintf(fid,'year, month, day, m3/s\n');
fclose(fid);
dlmwrite('Kongsbreen_discharge.csv',Q_kongsbreen,'-append');

%csvwrite('Kongsvegen_discharge.csv',{'year','month','day','m3/s'});
%csvwrite('Kongsvegen_discharge.csv',Q_kongsvegen,1,0);
%csvwrite('Kronebreen_discharge.csv',Q_kronebreen);
%csvwrite('Kongsbreen_discharge.csv',Q_kongsbreen);

% Get grid dimensions
xi = length(ncread(gname,'xi_rho'));
eta = length(ncread(gname,'eta_rho'));

% Get grid
lat = ncread(gname,'lat_rho');
lon = ncread(gname,'lon_rho');
theta = ncread(gname,'angle');
mask_rho = ncread(gname,'mask_rho');
depth = ncread(gname,'h');

% Kongsvegen
i1 = B(347,2) +1;  j1 = B(347,3) +1;  depth(i1,j1)
% Kronebreen
i3 = B(351,2) +1;  j3 = B(351,3) +1;  depth(i3,j3)
i4 = 81; j4 = 73;  depth(i4,j4)
% Kongsbreen
i5 = B(362,2) +1;  j5 = B(362,3) +1; depth(i5,j5)
i6 = 147; j6 = 71; depth(i6,j6)

first = 1;
for fn=1:nfiles
  % Read from model
  %rfil = [rdir,prename,num2str(fn,'%04d'),'.nc'];
  if t0+fn <= datenum([2007 06 1 0 0 0])
    rfil = [rdir,prename,datestr(t0+fn,'yyyymmdd'),'12'];
  else
      rfil = [rdir,prename,datestr(t0+fn,'yyyymmdd'),'11'];
  end
  
  if first
    % Find vertical depths
    Vtransform = ncread(rfil,'Vtransform');
    Vstretching = ncread(rfil,'Vstretching');
    theta_s = ncread(rfil,'theta_s');
    theta_b = ncread(rfil,'theta_b');
    hc = ncread(rfil,'Tcline');
    %hc = 5.;  %ncread(rfil,'Tcline');
    h0 = ncread(rfil,'h');
    s_rho = ncread(rfil,'s_rho'); N = length(s_rho); clear s_rho
    [z0] = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 1, h0, 0.); z0 = squeeze(abs(z0));
    [z0_w] = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 5, h0, 0.); z0_w = squeeze(abs(z0_w));    
    % clear Vtransform Vstretching theta_s theta_b hc
  end  % if first

  tid = ncread(rfil,'ocean_time'); tid = tid/(24*3600); dato = gregorian(tid + julian(1948,1,1));
  
  temp = ncread(rfil,'temp');
  salt = ncread(rfil,'salt');
  TEMP_1(fn,:) = temp(i4,j4,:);
  SALT_1(fn,:) = salt(i4,j4,:);  
  TEMP_2(fn,:) = temp(i6,j6,:);
  SALT_2(fn,:) = salt(i6,j6,:);
  TEMP_KoB(fn,:) = temp(i5,j5,:);
  SALT_KoB(fn,:) = salt(i5,j5,:);  
  TEMP_KrB(fn,:) = temp(i3,j3,:);
  SALT_KrB(fn,:) = salt(i3,j3,:);  
  TEMP_KoV(fn,:) = temp(i1,j1,:);
  SALT_KoV(fn,:) = salt(i1,j1,:);  
end

z0_Kongsvegen = [squeeze(z0(i1,j1,:)); NaN];
z0_w_Kongsvegen = squeeze(z0_w(i1,j1,:));
Kongsvegen_depth = [z0_w_Kongsvegen z0_Kongsvegen];
fid = fopen('Kongsvegen_depth.csv','w');
fprintf(fid,'H_grid, H_cell_center\n');
fclose(fid);
dlmwrite('Kongsvegen_depth.csv',Kongsvegen_depth,'-append');

z0_Kronebreen = [squeeze(z0(i3,j3,:)); NaN];
z0_w_Kronebreen = squeeze(z0_w(i3,j3,:));
Kronebreen_depth = [z0_w_Kronebreen z0_Kronebreen];
fid = fopen('Kronebreen_depth.csv','w');
fprintf(fid,'H_grid, H_cell_center\n');
fclose(fid);
dlmwrite('Kronebreen_depth.csv',Kronebreen_depth,'-append');

z0_Kronebreen_Ambient = [squeeze(z0(i4,j4,:)); NaN];
z0_w_Kronebreen_Ambient = squeeze(z0_w(i4,j4,:));
Kronebreen_depth_ambient = [z0_w_Kronebreen_Ambient z0_Kronebreen_Ambient];
fid = fopen('Kronebreen_depth_ambient.csv','w');
fprintf(fid,'H_grid, H_cell_center\n');
fclose(fid);
dlmwrite('Kronebreen_depth_ambient.csv',Kronebreen_depth_ambient,'-append');

z0_Kongsbreen = [squeeze(z0(i5,j5,:)); NaN];
z0_w_Kongsbreen = squeeze(z0_w(i5,j5,:));
Kongsbreen_depth = [z0_w_Kongsbreen z0_Kongsbreen];
fid = fopen('Kongsbreen_depth.csv','w');
fprintf(fid,'H_grid, H_cell_center\n');
fclose(fid);
dlmwrite('Kongsbreen_depth.csv',Kongsbreen_depth,'-append');

z0_Kongsbreen_Ambient = [squeeze(z0(i6,j6,:)); NaN];
z0_w_Kongsbreen_Ambient = squeeze(z0_w(i6,j6,:));
Kongsbreen_depth_ambient = [z0_w_Kongsbreen_Ambient z0_Kongsbreen_Ambient];
fid = fopen('Kongsbreen_depth_ambient.csv','w');
fprintf(fid,'H_grid, H_cell_center\n');
fclose(fid);
dlmwrite('Kongsbreen_depth_ambient.csv',Kongsbreen_depth_ambient,'-append');

%TEMP_Kronebreen = [z0_Kronebreen_Ambient TEMP_1'];
%SALT_Kronebreen = [z0_Kronebreen_Ambient SALT_1'];
%TEMP_Kongsbreen = [z0_Kongsbreen_Ambient TEMP_2'];
%SALT_Kongsbreen = [z0_Kongsbreen_Ambient SALT_2'];

dlmwrite('Kronebreen_Ambient_Temperature.csv',TEMP_1');
dlmwrite('Kronebreen_Ambient_Salinity.csv',SALT_1');
dlmwrite('Kongsbreen_Ambient_Temperature.csv',TEMP_2');
dlmwrite('Kongsbreen_Ambient_Salinity.csv',SALT_2');

%--------------------
% Convert to Absolute Salinity and Conservative Temperature, 
% calculate hydrostatic pressure and density.

% Get pressure from water depth (negative in ocean)
%P = gsw_p_from_z( -z0, lat(1) );

% Absolute Saliity
%SA = gsw_SA_from_SP( salt_mean, P, lon(1), lat(1) );

% Conservative Temperature
%CT = gsw_CT_from_pt( SA, temp_mean );

% Density
%RHO = gsw_rho( SA, CT, P );

%--------------------
lat1 = lat(i1,j1); lon1 = lon(i1,j1);
%
z_KoV = -squeeze(z0(i1,j1,:));
P_KoV = gsw_p_from_z( z_KoV, lat1 );
SA_KoV = gsw_SA_from_SP( SALT_KoV, P_KoV, lon1, lat1 );
CT_KoV = gsw_CT_from_pt( SA_KoV, TEMP_KoV );
RHO_KoV = gsw_rho( SA_KoV, CT_KoV, P_KoV );
save('Kongsvegen_profile.mat', 'z_KoV', 'P_KoV', 'SA_KoV', 'CT_KoV', 'RHO_KoV');
%
z_KrB = -squeeze(z0(i3,j3,:));
P_KrB = gsw_p_from_z( z_KrB, lat1 );
SA_KrB = gsw_SA_from_SP( SALT_KrB, P_KrB, lon1, lat1 );
CT_KrB = gsw_CT_from_pt( SA_KrB, TEMP_KrB );
RHO_KrB = gsw_rho( SA_KrB, CT_KrB, P_KrB );
save('Kronebreen_profile.mat', 'z_KrB', 'P_KrB', 'SA_KrB', 'CT_KrB', 'RHO_KrB');
%
z_KrB_amb = -squeeze(z0(i4,j4,:));
P_KrB_amb = gsw_p_from_z( z_KrB_amb, lat1 );
SA_KrB_amb = gsw_SA_from_SP( SALT_1, P_KrB_amb, lon1, lat1 );
CT_KrB_amb = gsw_CT_from_pt( SA_KrB_amb, TEMP_1 );
RHO_KrB_amb = gsw_rho( SA_KrB_amb, CT_KrB_amb, P_KrB_amb );
save('Kronebreen_amb_profile.mat', 'z_KrB_amb', 'P_KrB_amb', 'SA_KrB_amb', 'CT_KrB_amb', 'RHO_KrB_amb');
%
z_KoB = -squeeze(z0(i5,j5,:));
P_KoB = gsw_p_from_z( z_KoB, lat1 );
SA_KoB = gsw_SA_from_SP( SALT_KoB, P_KoB, lon1, lat1 );
CT_KoB = gsw_CT_from_pt( SA_KoB, TEMP_KoB );
RHO_KoB = gsw_rho( SA_KoB, CT_KoB, P_KoB );
save('Kongsbreen_profile.mat', 'z_KoB', 'P_KoB', 'SA_KoB', 'CT_KoB', 'RHO_KoB');
%
z_KoB_amb = -squeeze(z0(i6,j6,:));
P_KoB_amb = gsw_p_from_z( z_KoB_amb, lat1 );
SA_KoB_amb = gsw_SA_from_SP( SALT_2, P_KoB_amb, lon1, lat1 );
CT_KoB_amb = gsw_CT_from_pt( SA_KoB_amb, TEMP_2 );
RHO_KoB_amb = gsw_rho( SA_KoB_amb, CT_KoB_amb, P_KoB_amb );
save('Kongsbreen_amb_profile.mat', 'z_KoB_amb', 'P_KoB_amb', 'SA_KoB_amb', 'CT_KoB_amb', 'RHO_KoB_amb');
%