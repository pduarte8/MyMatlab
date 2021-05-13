clear all; close all;

% Define experiment
for expno=1:2

switch expno
case 1; exp = 'surface';
case 2; exp = 'subglacial';
end

% Interpolate all values to these depths
zout = [1,5,10:10:80];  %,150,200,250];

% Define borders
%lonmin = 11.+10/60; lonmax = 13.; latmin = 78.+50/60; latmax = 79.+25/60;  % ALL
lonmin = 12.+10/60; lonmax = 12.+40/60; latmin = 78.+51/60; latmax = 79.+1/60;  % Inner

% Read obs details
info = load('overview.txt');
lon = info(:,2); lat = info(:,3);
yy = info(:,4); mm = info(:,5); dd = info(:,6); hh = info(:,7);
clear info

% Read, sort and store obs (assume that downcast is stored on files, only)
for i=1:length(lon)
  obs = load(['Data/ctd_',num2str(i,'%04d'),'.asc']);
  zlev_unsorted = obs(:,1);
  [zlev_sorted,IA,IC] = unique(zlev_unsorted);
  nvalues(i) = length(zlev_sorted);
  zlev(i,1:nvalues(i)) = zlev_sorted;
  osalt(i,1:nvalues(i)) = obs(IA,2);
  otemp(i,1:nvalues(i)) = obs(IA,3);
  clear obs IA IC zlev_unsorted zlev_sorted
end
clear obs i

% Read model land mask, used for plotting
gfil = '/cluster/work/users/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/kongsfjorden_160m_grid_present.nc';
mask_rho = ncread(gfil,'mask_rho');
lon_rho = ncread(gfil,'lon_rho');
lat_rho = ncread(gfil,'lat_rho');
clear gfil

% Read model values
inp = load(['ctd_stations_Kongsfjorden_',exp,'.asc']);
YY = inp(:,1); MM = inp(:,2); DD = inp(:,3); HH = inp(:,4); LON = inp(:,5); LAT = inp(:,6);
ZLEV = inp(:,7:41); MSALT = inp(:,42:76); MTEMP = inp(:,77:111);
clear inp

% Model input defines which CTD-casts are inside fjord, interpolate to observational depths and exclude obs outside fjord
YYYYMMDDHH = YY*1000000 + MM*10000 + DD*100 + HH;
yyyymmddhh = yy*1000000 + mm*10000 + dd*100 + hh;

for p=1:length(YY)
  n = find(YYYYMMDDHH(p) == yyyymmddhh & LON(p) == lon & LAT(p) == lat);
  ms = MSALT(p,:); mt = MTEMP(p,:); mz = ZLEV(p,:);
  os = osalt(n,1:nvalues(n)); ot = otemp(n,1:nvalues(n)); oz = zlev(n,1:nvalues(n));
  clear n
  obssalt(p,1:length(zout)) = interp1(oz,os,zout);
  obstemp(p,1:length(zout)) = interp1(oz,ot,zout);
  modsalt(p,1:length(zout)) = interp1(mz,ms,zout);
  modtemp(p,1:length(zout)) = interp1(mz,mt,zout);
  dato(p,1) = YY(p); dato(p,2) = MM(p); dato(p,3) = DD(p); dato(p,4) = HH(p);
  clear ms mt mz os ot oz
end
clear YYYYMMDDHH YY MM DD HH yyyymmddhh yy mm dd hh p MSALT MTEMP ZLEV osalt otemp zlev lon lat nvalues

% Exlude obs outside chosen focus area
n = find(LON >= lonmin & LON <= lonmax & LAT >= latmin & LAT <= latmax);
LON = LON(n); LAT = LAT(n); dato = dato(n,:); obssalt = obssalt(n,:); modsalt = modsalt(n,:);
                                              obstemp = obstemp(n,:); modtemp = modtemp(n,:);
clear n

% Make error statistics
for k=1:length(zout)
  sbias(k,expno) = nanmean(modsalt(:,k) - obssalt(:,k));
  tbias(k,expno) = nanmean(modtemp(:,k) - obstemp(:,k));
end

OBSSALT(:,:,expno) = obssalt; OBSTEMP(:,:,expno) = obstemp;
MODSALT(:,:,expno) = modsalt; MODTEMP(:,:,expno) = modtemp;

clear dato k LAT lat_rho latmax latmin LON lon_rho lonmax lonmin mask_rho obssalt obstemp modsalt modtemp

end  % expno

% Plot model error with depth
fsize = 30;  % Fontsize
figure('position',[5 5 70 150])
clf
plot(sbias(:,1),-zout,'b-','LineWidth',4); hold on;
plot(sbias(:,2),-zout,'r-','LineWidth',4); hold on;
legend('PSD','PTG','Location','SouthWest');
for k=1:length(zout)
  errors = MODSALT(:,k,1)-OBSSALT(:,k,1);                        % All model deviations from obs
  n = find(~isnan(errors)); errors = errors(n);                  % Remove NaNs
  err_neg = prctile(errors,10.); err_pos = prctile(errors,90.);  % Lower and upper value of error bar
  vals = [err_neg:0.01:err_pos];                                 % Make vector denoting the error bar
  plot(vals,repmat(-zout(k),length(vals)),'b-','LineWidth',2.0); hold on;
  plot([err_neg,err_pos],[-zout(k),-zout(k)],'b+','MarkerSize',5); hold on;
  clear errors n err_neg err_pos vals
end
for k=1:length(zout)
  errors = MODSALT(:,k,2)-OBSSALT(:,k,2);                        % All model deviations from obs
  n = find(~isnan(errors)); errors = errors(n);                  % Remove NaNs
  err_neg = prctile(errors,10.); err_pos = prctile(errors,90.);  % Lower and upper value of error bar
  vals = [err_neg:0.01:err_pos];                                 % Make vector denoting the error bar
  plot(vals,repmat(-zout(k),length(vals)),'r-','LineWidth',2.0); hold on;
  plot([err_neg,err_pos],[-zout(k),-zout(k)],'r+','MarkerSize',7); hold on;
  clear errors n err_neg err_pos vals
end
plot(zout*0,-zout,'k-','LineWidth',1.5); hold on;
set(gca,'FontSize',fsize,'FontWeight','bold');
xlabel('Salinity error (mod-obs)','FontSize',fsize,'FontWeight','bold');
ylabel('Depth (m)','FontSize',fsize,'FontWeight','bold');
set(gcf,'PaperPositionMode','auto');
print('-dpng','-r300',['kongsfjorden_errorswithdepth_CTD_salt.png']);

figure('position',[5 5 70 150])
clf
plot(tbias(:,1),-zout,'b-','LineWidth',4); hold on;
plot(tbias(:,2),-zout,'r-','LineWidth',4); hold on;
legend('PSD','PTG','Location','SouthWest');
for k=1:length(zout)
  errors = MODTEMP(:,k,1)-OBSTEMP(:,k,1);                        % All model deviations from obs
  n = find(~isnan(errors)); errors = errors(n);                  % Remove NaNs
  err_neg = prctile(errors,10.); err_pos = prctile(errors,90.);  % Lower and upper value of error bar
  vals = [err_neg:0.01:err_pos];                                 % Make vector denoting the error bar
  plot(vals,repmat(-zout(k),length(vals)),'b-','LineWidth',2.0); hold on;
  plot([err_neg,err_pos],[-zout(k),-zout(k)],'b+','MarkerSize',5); hold on;
  clear errors n err_neg err_pos vals
end
for k=1:length(zout)
  errors = MODTEMP(:,k,2)-OBSTEMP(:,k,2);                        % All model deviations from obs
  n = find(~isnan(errors)); errors = errors(n);                  % Remove NaNs
  err_neg = prctile(errors,10.); err_pos = prctile(errors,90.);  % Lower and upper value of error bar
  vals = [err_neg:0.01:err_pos];                                 % Make vector denoting the error bar
  plot(vals,repmat(-zout(k),length(vals)),'r-','LineWidth',2.0); hold on;
  plot([err_neg,err_pos],[-zout(k),-zout(k)],'r+','MarkerSize',7); hold on;
  clear errors n err_neg err_pos vals
end
plot(zout*0,-zout,'k-','LineWidth',1.5); hold on;
set(gca,'FontSize',fsize,'FontWeight','bold');
xlabel('Temperature error (mod-obs)','FontSize',fsize,'FontWeight','bold');
ylabel('Depth (m)','FontSize',fsize,'FontWeight','bold');
set(gcf,'PaperPositionMode','auto');
print('-dpng','-r300',['kongsfjorden_errorswithdepth_CTD_temp.png']);
