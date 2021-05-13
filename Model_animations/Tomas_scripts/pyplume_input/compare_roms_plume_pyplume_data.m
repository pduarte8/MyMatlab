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

rdir = 'D:\Data\roms_run\Kongsfjorden-160m\K160_v2\present_plume\daily_avg\';
prename = 'kongsfjorden_160m_avg.nc4_';
%nfiles = 10;
%nfiles = 90;
nfiles = 121;
t0 = datenum([2007 05 31 0 0 0]);  % Last day before beginning of data record

TEMP_1 = zeros(nfiles,35);
SALT_1 = zeros(nfiles,35);
TEMP_2 = zeros(nfiles,35);
SALT_2 = zeros(nfiles,35);
TEMP_KoB = zeros(nfiles,35);
SALT_KoB = zeros(nfiles,35);
TEMP_KrB = zeros(nfiles,35);
SALT_KrB = zeros(nfiles,35);
TEMP_KoV = zeros(nfiles,35);
SALT_KoV = zeros(nfiles,35);

u1_KoB = zeros(nfiles,35);
u2_KoB = zeros(nfiles,35);
v1_KoB = zeros(nfiles,35);
v2_KoB = zeros(nfiles,35);
vf_u_KoB = zeros(nfiles,35);
vf_v_KoB = zeros(nfiles,35);
vf_KoB = zeros(nfiles,35);

u1_KrB = zeros(nfiles,35);
u2_KrB = zeros(nfiles,35);
v1_KrB = zeros(nfiles,35);
v2_KrB = zeros(nfiles,35);
vf_u_KrB = zeros(nfiles,35);
vf_v_KrB = zeros(nfiles,35);
vf_KrB = zeros(nfiles,35);

u1_KoV = zeros(nfiles,35);
u2_KoV = zeros(nfiles,35);
v1_KoV = zeros(nfiles,35);
v2_KoV = zeros(nfiles,35);
vf_u_KoV = zeros(nfiles,35);
vf_v_KoV = zeros(nfiles,35);
vf_KoV = zeros(nfiles,35);

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
Q_kongsvegen = squeeze(Q_flux(i1:i2,347));
Q_kronebreen = squeeze(Q_flux(i1:i2,351));
Q_kongsbreen = squeeze(Q_flux(i1:i2,362));

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
i_KoV = B(347,2) +1;  j_KoV = B(347,3) +1;  %depth(i_KoV,j_KoV)
% Kronebreen
i_KrB = B(351,2) +1;  j_KrB = B(351,3) +1;  %depth(i_KrB,j_KrB)
i4 = 81; j4 = 73;  %depth(i4,j4)
% Kongsbreen
i_KoB = B(362,2) +1;  j_KoB = B(362,3) +1; %depth(i_KoB,j_KoB)
i6 = 147; j6 = 71; %depth(i6,j6)

first = 1;
for fn=1:nfiles
  % Read from model
  %rfil = [rdir,prename,num2str(fn,'%04d'),'.nc'];
  %if t0+fn <= datenum([2007 06 1 0 0 0])
    rfil = [rdir,prename,datestr(t0+fn,'yyyymmdd'),'12'];
  %else
  %    rfil = [rdir,prename,datestr(t0+fn,'yyyymmdd'),'11'];
  %end
  disp(rfil);
  
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

    z0_Kongsvegen = squeeze(z0(i_KoV,j_KoV,:));
    z0_w_Kongsvegen = squeeze(z0_w(i_KoV,j_KoV,:));

    z0_Kronebreen = squeeze(z0(i_KrB,j_KrB,:));
    z0_w_Kronebreen = squeeze(z0_w(i_KrB,j_KrB,:));

    z0_Kronebreen_Ambient = squeeze(z0(i4,j4,:));
    z0_w_Kronebreen_Ambient = squeeze(z0_w(i4,j4,:));

    z0_Kongsbreen = squeeze(z0(i_KoB,j_KoB,:));
    z0_w_Kongsbreen = squeeze(z0_w(i_KoB,j_KoB,:));

    z0_Kongsbreen_Ambient = squeeze(z0(i6,j6,:));
    z0_w_Kongsbreen_Ambient = squeeze(z0_w(i6,j6,:));

    first = 0;
  end  % if first

  tid = ncread(rfil,'ocean_time'); tid = tid/(24*3600); dato = gregorian(tid + julian(1948,1,1));
  
  temp = ncread(rfil,'temp');
  salt = ncread(rfil,'salt');
  TEMP_1(fn,:) = temp(i4,j4,:);
  SALT_1(fn,:) = salt(i4,j4,:);  
  TEMP_2(fn,:) = temp(i6,j6,:);
  SALT_2(fn,:) = salt(i6,j6,:);
  TEMP_KoV(fn,:) = temp(i_KoV,j_KoV,:);
  SALT_KoV(fn,:) = salt(i_KoV,j_KoV,:);
  TEMP_KrB(fn,:) = temp(i_KrB,j_KrB,:);
  SALT_KrB(fn,:) = salt(i_KrB,j_KrB,:);
  TEMP_KoB(fn,:) = temp(i_KoB,j_KoB,:);
  SALT_KoB(fn,:) = salt(i_KoB,j_KoB,:);
  
  u = ncread(rfil,'u');
  v = ncread(rfil,'v');
%   u1_KoB(fn,:) = u(i_KoB,j_KoB,:);
%   u2_KoB(fn,:) = u(i_KoB+1,j_KoB,:);
%   v1_KoB(fn,:) = v(i_KoB,j_KoB,:);
%   v2_KoB(fn,:) = v(i_KoB,j_KoB+1,:);
% 
%   u1_KrB(fn,:) = u(i_KrB,j_KrB,:);
%   u2_KrB(fn,:) = u(i_KrB+1,j_KrB,:);
%   v1_KrB(fn,:) = v(i_KrB,j_KrB,:);
%   v2_KrB(fn,:) = v(i_KrB,j_KrB+1,:);
% 
%   u1_KoV(fn,:) = u(i_KoV,j_KoV,:);
%   u2_KoV(fn,:) = u(i_KoV+1,j_KoV,:);
%   v1_KoV(fn,:) = v(i_KoV,j_KoV,:);
%   v2_KoV(fn,:) = v(i_KoV,j_KoV+1,:);
%     
%   for k=1:35
%       vf_u_KoB(fn,k) = u2_KoB(fn,k) - u1_KoB(fn,k);
%       vf_v_KoB(fn,k) = v2_KoB(fn,k) - v1_KoB(fn,k);
%       vf_KoB(fn,k) = (vf_u_KoB(fn,k) + vf_v_KoB(fn,k)) * (z0_w_Kongsbreen(k) - z0_w_Kongsbreen(k+1)) * 160.0;
% 
%       vf_u_KrB(fn,k) = u2_KrB(fn,k) - u1_KrB(fn,k);
%       vf_v_KrB(fn,k) = v2_KrB(fn,k) - v1_KrB(fn,k);
%       vf_KrB(fn,k) = (vf_u_KrB(fn,k) + vf_v_KrB(fn,k)) * (z0_w_Kronebreen(k) - z0_w_Kronebreen(k+1)) * 160.0;
% 
%       vf_u_KoV(fn,k) = u2_KoV(fn,k) - u1_KoV(fn,k);
%       %vf_v_KoV(fn,k) = v2_KoV(fn,k) - v1_KoV(fn,k);
%       %vf_KoV(fn,k) = (vf_u_KoV(fn,k) + vf_v_KoV(fn,k)) * (z0_w_Kongsvegen(k) - z0_w_Kongsvegen(k+1)) * 160.0;
%       vf_KoV(fn,k) = (vf_u_KoV(fn,k) - v1_KoV(fn,k)) * (z0_w_Kongsvegen(k) - z0_w_Kongsvegen(k+1)) * 160.0;
%   end
% end

  u1_KoB(fn,:) = u(i_KoB-1,j_KoB,:);   u1_KoB(isnan(u1_KoB))=0;
  u2_KoB(fn,:) = u(i_KoB,j_KoB,:);     u2_KoB(isnan(u2_KoB))=0;
  v1_KoB(fn,:) = v(i_KoB,j_KoB-1,:);   v1_KoB(isnan(v1_KoB))=0;
  v2_KoB(fn,:) = v(i_KoB,j_KoB,:);     v2_KoB(isnan(v2_KoB))=0;

  u1_KrB(fn,:) = u(i_KrB-1,j_KrB,:);   u1_KrB(isnan(u1_KrB))=0;
  u2_KrB(fn,:) = u(i_KrB,j_KrB,:);     u2_KrB(isnan(u2_KrB))=0;
  v1_KrB(fn,:) = v(i_KrB,j_KrB-1,:);   v1_KrB(isnan(v1_KrB))=0;
  v2_KrB(fn,:) = v(i_KrB,j_KrB,:);     v2_KrB(isnan(v2_KrB))=0;

  u1_KoV(fn,:) = u(i_KoV-1,j_KoV,:);   u1_KoV(isnan(u1_KoV))=0;
  u2_KoV(fn,:) = u(i_KoV,j_KoV,:);     u2_KoV(isnan(u2_KoV))=0;
  v1_KoV(fn,:) = v(i_KoV,j_KoV-1,:);   v1_KoV(isnan(v1_KoV))=0;
  v2_KoV(fn,:) = v(i_KoV,j_KoV,:);     v2_KoV(isnan(v2_KoV))=0;
    
  for k=1:35
      vf_u_KoB(fn,k) = u2_KoB(fn,k) - u1_KoB(fn,k);
      vf_v_KoB(fn,k) = v2_KoB(fn,k) - v1_KoB(fn,k);
      vf_KoB(fn,k) = (vf_u_KoB(fn,k) + vf_v_KoB(fn,k)) * (z0_w_Kongsbreen(k) - z0_w_Kongsbreen(k+1)) * 160.0;

      vf_u_KrB(fn,k) = u2_KrB(fn,k) - u1_KrB(fn,k);
      vf_v_KrB(fn,k) = v2_KrB(fn,k) - v1_KrB(fn,k);
      vf_KrB(fn,k) = (vf_u_KrB(fn,k) + vf_v_KrB(fn,k)) * (z0_w_Kronebreen(k) - z0_w_Kronebreen(k+1)) * 160.0;

      vf_u_KoV(fn,k) = u2_KoV(fn,k) - u1_KoV(fn,k);
      vf_v_KoV(fn,k) = v2_KoV(fn,k) - v1_KoV(fn,k);
      vf_KoV(fn,k) = (vf_u_KoV(fn,k) + vf_v_KoV(fn,k)) * (z0_w_Kongsvegen(k) - z0_w_Kongsvegen(k+1)) * 160.0;
  end
end

% Load plume data
Zpp = csvread('Kongsbreen_output_86m\kongsbreen_results_z.csv');
Tpp = csvread('Kongsbreen_output_86m\kongsbreen_results_t_p.csv');
Spp = csvread('Kongsbreen_output_86m\kongsbreen_results_s_p.csv');

% for i=1:122
%     plot(Spp(i,:),Zpp(i,:)-Zpp(i,end),'b-',SALT_KoB(i,:),-z0_Kongsbreen,'r-',SALT_1(i,:),-z0_Kongsbreen_Ambient,'g-')
%     pause(0.25)
% end

%% Use only inflow
% vf_KoB_i = vf_KoB;
% vf_KoB_o = vf_KoB;
% vf_KrB_i = vf_KrB;
% vf_KrB_o = vf_KrB;
% vf_KoV_i = vf_KoV;
% vf_KoV_o = vf_KoV;
% 
% vf_KoB_i(vf_KoB_i < 0) = 0;
% vf_KoB_o(vf_KoB_o > 0) = 0;
% vf_KrB_i(vf_KrB_i < 0) = 0;
% vf_KrB_o(vf_KrB_o > 0) = 0;
% vf_KoV_i(vf_KoV_i < 0) = 0;
% vf_KoV_o(vf_KoV_o > 0) = 0;
%
vf_KoB(vf_KoB < 0) = 0;
vf_KrB(vf_KrB < 0) = 0;
vf_KoV(vf_KoV < 0) = 0;
%

% vf_KoB_i_sum = sum(vf_KoB_i(:,:),2);
% vf_KoB_o_sum = -sum(vf_KoB_o(:,:),2);
% vf_KrB_i_sum = sum(vf_KrB_i(:,:),2);
% vf_KrB_o_sum = -sum(vf_KrB_o(:,:),2);
% vf_KoV_i_sum = sum(vf_KoV_i(:,:),2);
% vf_KoV_o_sum = -sum(vf_KoV_o(:,:),2);
vf_KoB_sum = sum(vf_KoB(:,:),2);
vf_KrB_sum = sum(vf_KrB(:,:),2);
vf_KoV_sum = sum(vf_KoV(:,:),2);

DisEnt_KoB = csvread('Kongsbreen_output_86m\kongsbreen_results_disch_ent_n1.csv',1,0);
figure
scatter(DisEnt_KoB(:,5),DisEnt_KoB(:,9),5.0,'o');
hold on;
scatter(DisEnt_KoB(1:121,5),vf_KoB_sum,5.0,'o');
%scatter(DisEnt_KoB(1:121,5),vf_KoB_o_sum,5.0,'o');
legend('pyplume','ROMS','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsbreen (plume)');

DisEnt_KoB = csvread('Kongsbreen_output_86m\kongsbreen_results_disch_ent_n1.csv',1,0);
figure
scatter(DisEnt_KoB(:,5),DisEnt_KoB(:,9),5.0,'o');
hold on;
scatter(DisEnt_KoB(1:121,5),vf_KoB_sum/10,5.0,'o');
%scatter(DisEnt_KoB(1:121,5),vf_KoB_o_sum,5.0,'o');
legend('pyplume','ROMS iceplume / 10','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsbreen (plume)');

DisEnt_KrB = csvread('Kronebreen_output\kronebreen_results_disch_ent_nb.csv',1,0);
figure
scatter(DisEnt_KrB(:,5),DisEnt_KrB(:,9),5.0,'o');
hold on;
scatter(DisEnt_KrB(1:121,5),vf_KrB_sum,5.0,'o');
%scatter(DisEnt_KrB(1:121,5),vf_KrB_o_sum,5.0,'o');
legend('pyplume','ROMS','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kronebreen (plume)');

DisEnt_KoV = csvread('Kongsvegen_output\kongsvegen_results_disch_ent_nb.csv',1,0);
figure
scatter(DisEnt_KoV(:,5),DisEnt_KoV(:,9),5.0,'o');
hold on;
scatter(DisEnt_KoV(1:121,5),vf_KoV_sum,5.0,'o');
%scatter(DisEnt_KoV(1:121,5),vf_KoV_o_sum,5.0,'o');
legend('pyplume','ROMS','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsvegen (plume)');

vf_KoB_plm = vf_KoB_sum;
vf_KrB_plm = vf_KrB_sum;
vf_KoV_plm = vf_KoV_sum;

save('VolumeFlux_plm.mat', 'vf_KoB_plm', 'vf_KrB_plm', 'vf_KoV_plm');