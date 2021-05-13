clear all; close all;
addpath('/cluster/home/pduarte/matlabmfiles/ROMS/matlab/tools');
%addpath('/cluster/home/pduarte/matlabmfiles/ROMS/matlab/mexcdf5/');
addpath('/cluster/home/pduarte/matlabmfiles/ROMS/matlab/mexcdf/mexnc/');
% Find S,T-values from model archive based on a list of lon,lat,z-positions

exp = 'present_subglacial';
%nfiles = 121;
% Define subdomain to reduce memory use
i0 = 1; i1 = 400; j0 = 1; j1 = 600;

% Model input - static fields
rdir = ['/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/ThirdRoundOfSimulations/Sim12_repeated';
%rdir = ['/cluster/home/pduarte/matlabmfiles/Temp/'];
rdir = ['/cluster/shared/arcticfjord/run_Kongsfjorden-160m_present_subglacial/temp/'];
prename = 'ocean_his.nc_';
rfil = [rdir,prename,'2007060100-2007060200']
%lon_rho = ncread(rfil,'lon_rho',[i0 j0],[i1-i0+1 j1-j0+1]);
lon = ncread(rfil,'lon_rho');
lon(1,1)
lon_rho = lon(i0:i1-i0+1, j0:j1-j0+1);
%lat_rho = ncread(rfil,'lat_rho',[i0 j0],[i1-i0+1 j1-j0+1]);
lat = ncread(rfil,'lat_rho');
lat_rho = lat(i0:i1-i0+1,j0:j1-j0+1); 
%Vtransform = ncread(rfil,'Vtransform'); Vstretching = ncread(rfil,'Vstretching');
%theta_s = ncread(rfil,'theta_s'); theta_b = ncread(rfil,'theta_b'); hc = ncread(rfil,'hc');
Vtransform = 2.0; Vstretching = 2.0;
theta_s = 8.0; theta_b = 0.1; hc = 20.0;
%h = ncread(rfil,'h',[i0 j0],[i1-i0+1 j1-j0+1]); 
hh = ncread(rfil,'h');
h = hh(i0:i1-i0+1,j0:j1-j0+1);
s_rho = ncread(rfil,'s_rho'); N = length(s_rho);
clear s_rho rfil

% Read positions
pfil = 'overview.txt';
inp = load(pfil);
yy = inp(:,4); mm = inp(:,5); dd = inp(:,6); hh = inp(:,7); lat = inp(:,3); lon = inp(:,2);
clear inp pfil

%ofil = ['/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc_shared/Sim5_PTG0/Temp/',exp_job,'.asc'];
ofil = ['/cluster/home/pduarte/matlabmfiles/Temp/',exp,'.asc'];

% Print to file
fid = fopen(ofil,'w');

% Loop over positions
for p=1:length(lon)

  % Find correct file, i.e. date
  t0 = yy(p)*1000000 + mm(p)*10000 + dd(p)*100 + 1;
  t1 = yy(p)*1000000 + mm(p)*10000 + (dd(p)+1)*100;
  %found = false
  %i = 0;
  %while i < nfiles & found == false
  %    i = i + 1;
  %    info = ncinfo(
  %    fname = n
  
  rfil = [rdir,prename,int2str(t0),'-',int2str(t1)];
  tid = ncread(rfil,'ocean_time');
  dato = gregorian(tid/(24*3600) + julian(1948,1,1,0));
  tstep = find(hh(p) == dato(:,4));
  clear tid dato
  tstep
  %salt3d = ncread(rfil,'salt',[i0 j0 1 tstep],[i1-i0+1 j1-j0+1 Inf 1]);
  %temp3d = ncread(rfil,'temp',[i0 j0 1 tstep],[i1-i0+1 j1-j0+1 Inf 1]);
  salt3dd = ncread(rfil,'salt');
  temp3dd = ncread(rfil,'temp');
  salt3d = salt3dd(i0:i1-i0+1,j0:j1-j0+1,:,tstep);
  temp3d = temp3dd(i0:i1-i0+1,j0:j1-j0+1,:,tstep);
  clear tstep t0 t1

  % Interpolate horizontally
  for k=1:N
    F = scatteredInterpolant(reshape(lon_rho,numel(lon_rho),1),reshape(lat_rho,numel(lat_rho),1),reshape(salt3d(:,:,k),numel(lon_rho),1),'linear','none');
    saltp(k) = F(lon(p),lat(p));
    clear F
    F = scatteredInterpolant(reshape(lon_rho,numel(lon_rho),1),reshape(lat_rho,numel(lat_rho),1),reshape(temp3d(:,:,k),numel(lon_rho),1),'linear','none');
    tempp(k) = F(lon(p),lat(p));
    clear F
  end
  clear salt3d temp3d k

  if ~isnan(saltp(end))
    F = scatteredInterpolant(reshape(lon_rho,numel(lon_rho),1),reshape(lat_rho,numel(lat_rho),1),reshape(h,numel(lon_rho),1),'linear','none');
    totdepth = F(lon(p),lat(p));
    clear F
    [z] = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,1,totdepth,0.); z = squeeze(abs(z));
    % Print profile to file
    fprintf(fid,'%04d %02d %02d %02d %9.6f %9.6f', yy(p), mm(p), dd(p), hh(p), lon(p), lat(p));
    for k=length(z):-1:1; fprintf(fid,'%7.2f', z(k)); end  % k
    for k=length(z):-1:1; fprintf(fid,'%7.2f', saltp(k)); end  % k
    for k=length(z):-1:1; fprintf(fid,'%7.2f', tempp(k)); end  % k
    fprintf(fid,'\n');
    clear saltp tempp totdepth z
    disp(['Finished interpolating S and T to position no. ',int2str(p),' (',int2str(length(lon)),')']);
  else
    disp(['Position no. ',int2str(p),' (',int2str(length(lon)),') was outside model domain']);
  end

end  % p

fclose(fid);
