clear all; close all;
addpath /cluster/home/pduarte/matlabmfiles/m_map
% Define experiment
for expno=1:4

switch expno
case 1; exp = 'Sim12_present_subglacial.asc';
case 2; exp = 'Sim21_present_subglacial.asc';
case 3; exp = 'Sim23_present_subglacial.asc';    
case 4; exp = 'Sim25_present_subglacial.asc';    
end
% Interpolate all values to these depths
zout = [1,10:10:50];  %,150,200,250];

% Define borders
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
inp = load(exp);
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
LON = LON(n); LAT = LAT(n); dato = dato(n,:); obssalt = obssalt(n,:); modsalt = modsalt(n,:); obstemp = obstemp(n,:); modtemp = modtemp(n,:);
clear n

% Make error statistics
for k=1:length(zout)
  sbias(k) = nanmean(modsalt(:,k) - obssalt(:,k));
  tbias(k) = nanmean(modtemp(:,k) - obstemp(:,k));
end

% Compare arrays from obs and mod, if one is NaN, set the other to NaN as well
for s=1:size(obssalt,1)
  for k=1:length(zout)
    if isnan(obssalt(s,k)); modsalt(s,k) = NaN; end
    if isnan(modsalt(s,k)); obssalt(s,k) = NaN; end
    if isnan(obstemp(s,k)); modtemp(s,k) = NaN; end
    if isnan(modtemp(s,k)); obstemp(s,k) = NaN; end
  end
end

% Plot scatter and qq
for z=zout
  k = find(z == zout);
  min_S = floor(min(min(obssalt(:,k),modsalt(:,k)))); max_S = ceil(max(max(obssalt(:,k),modsalt(:,k))));
  min_T = floor(min(min(obstemp(:,k),modtemp(:,k)))); max_T = ceil(max(max(obstemp(:,k),modtemp(:,k))));
  x_S = [min_S:0.5:max_S];
  x_T = [min_T:1.0:max_T];
  [no_s x] = hist(obssalt(:,k),x_S); no_s = no_s/sum(no_s);
  [no_t x] = hist(obstemp(:,k),x_T); no_t = no_t/sum(no_t);
  [nm_s x] = hist(modsalt(:,k),x_S); nm_s = nm_s/sum(nm_s);
  [nm_t x] = hist(modtemp(:,k),x_T); nm_t = nm_t/sum(nm_t);

  figure(1)
  clf
  subplot(2,3,1)
  plot(obssalt(:,k),modsalt(:,k),'b*','MarkerSize',2); hold on;
  plot([min_S:0.1:max_S],[min_S:0.1:max_S],'k-'); hold on;
  axis([min_S max_S min_S max_S]); axis square;
  xlabel('Obs'); ylabel('Mod');
  text(min_S,max_S+0.1*(max_S-min_S),'Scatter','FontSize',10,'FontWeight','bold');
  text(min_S,min_S-0.5*(max_S-min_S),['All CTD casts at ',int2str(z),'m depth'],'FontSize',12,'FontWeight','bold');
  subplot(2,3,2)
  plot(sort(obssalt(:,k)),sort(modsalt(:,k)),'b*','MarkerSize',2); hold on;
  plot([min_S:0.1:max_S],[min_S:0.1:max_S],'k-'); hold on;
  axis([min_S max_S min_S max_S]); axis square;
  xlabel('Obs'); ylabel('Mod');
  text(min_S,max_S+0.1*(max_S-min_S),'QQ','FontSize',10,'FontWeight','bold');
  subplot(2,3,3)
  line(x_S,no_s,'Color','k','LineWidth',1.0,'LineStyle','-'); hold on;
  line(x_S,nm_s,'Color','b','LineWidth',1.0,'LineStyle','-'); hold on;
  axis square;
  xlabel('Salinity'); ylabel('Freq');
  legend(['Obs (',num2str(nanmean(obssalt(:,k)),'%.1f'),')'],['Mod (',num2str(nanmean(modsalt(:,k)),'%.1f'),')'],'Location','SouthOutside');
  subplot(2,3,4)
  plot(obstemp(:,k),modtemp(:,k),'b*','MarkerSize',2); hold on;
  plot([min_T:0.1:max_T],[min_T:0.1:max_T],'k-'); hold on;
  axis([min_T max_T min_T max_T]); axis square;
  xlabel('Obs'); ylabel('Mod');
  subplot(2,3,5)
  plot(sort(obstemp(:,k)),sort(modtemp(:,k)),'b*','MarkerSize',2); hold on;
  plot([min_T:0.1:max_T],[min_T:0.1:max_T],'k-'); hold on;
  axis([min_T max_T min_T max_T]); axis square;
  xlabel('Obs'); ylabel('Mod');
  subplot(2,3,6)
  line(x_T,no_t,'Color','k','LineWidth',1.0,'LineStyle','-'); hold on;
  line(x_T,nm_t,'Color','b','LineWidth',1.0,'LineStyle','-'); hold on;
  axis square;
  xlabel('Temperature (^oC)'); ylabel('Freq');
  legend(['Obs (',num2str(nanmean(obstemp(:,k)),'%.1f'),'^oC)'],['Mod (',num2str(nanmean(modtemp(:,k)),'%.1f'),'^oC)'],'Location','SouthOutside');
  print('-dpng','-r300',['kongsfjorden_statistics_CTD_',num2str(z,'%02d'),'m_',exp,'.png']);

  % Define colors on colorbar
  b = [32 96 255;32 191 255;42 255 255;255 255 84;255 255 84;255 191 0;255 138 0;255 0 0]/255;

  figure(2)
  clf
  colormap(b)
  subplot(1,2,1)
  m_proj('lambert','long',[lonmin-5/60 lonmax+5/60],'lat',[latmin latmax]);
  m_contour(lon_rho,lat_rho,mask_rho,[0.5 0.5],'k-'); hold on;
  m_scatter(LON,LAT,[30],modsalt(:,k)-obssalt(:,k),'filled'); hold on;
  caxis([-1.5 1.5]);
  cb = colorbar;
  cb.Label.String = 'Salinity';
  title(['Smod - Sobs (',int2str(z),'m)'],'FontSize',14,'FontWeight','bold');
  m_grid('box','fancy','tickdir','in');
  subplot(1,2,2)
  m_proj('lambert','long',[lonmin-5/60 lonmax+5/60],'lat',[latmin latmax]);
  m_contour(lon_rho,lat_rho,mask_rho,[0.5 0.5],'k-'); hold on;
  m_scatter(LON,LAT,[30],modtemp(:,k)-obstemp(:,k),'filled'); hold on;
  caxis([-1.5 1.5]);
  cb = colorbar;
  cb.Label.String = 'Temperature';
  title(['Tmod - Tobs (',int2str(z),'m)'],'FontSize',14,'FontWeight','bold');
  m_grid('box','fancy','tickdir','in');
  print('-dpng','-r300',['kongsfjorden_scattermap_CTD_',num2str(z,'%02d'),'m_',exp,'.png']);

  clear no_s nm_s no_t nm_t x
end

clearvars -except expno

end
