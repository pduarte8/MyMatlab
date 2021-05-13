clear;

% Plot parameters
sim_label = 'K160v2';
%
 case_label = 'subglacial';
% case_label = 'surface';
% case_label = 'future';
%
% month = 'June';
% month = 'July';
% month = 'Aug';
month = 'Sep';

% Depth range : water depths are negative
%               use large positive value e.g. 99 to include sea surface
 depth_range = [ -30   99];  range_label = '0to30m';
% depth_range = [ -30  -20];  range_label = '20to30m';
% depth_range = [ -55  -45];  range_label = '45to55m';
% depth_range = [ -80  -70];  range_label = '70to80m';
% depth_range = [-105  -95];  range_label = '95to105m';
% depth_range = [ -10   1000];  range_label = '0to10m'; %The one I used
% depth_range = [ -5   1000];  range_label = '0to5m'; 
% depth_range = [ -30   1000];  range_label = '0to30m'; 

label = sprintf('%s_%s_%s', case_label,month,range_label);

t0 = datenum([1948 1 1 0 0 0]);

N = 35;

% Grid file
switch case_label
    case 'subglacial'
        % Grid file
        gname = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/kongsfjorden_160m_grid_present.nc';
        % Data file directory
        dir_data = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/Sim13/';
        % file list
        file_list = sprintf('files_%s_%s_%s.txt', sim_label,case_label,month);
    case 'future'
        % Grid file
        gname = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_future/kongsfjorden_160m_grid_future.nc';
        % Data file directory
        dir_data = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_future/Sim6/';
        % file list
        file_list = sprintf('files_%s_%s_%s.txt', sim_label,case_label,month);
    case 'surface'
        % Grid file
        gname = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/kongsfjorden_160m_grid_present.nc';
        % Data file directory
        dir_data = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_present/Sim4/';
        % file list
        file_list = sprintf('files_%s_subglacial_%s.txt', sim_label,month);
    otherwise
        error('Not a valid case.')
end

prename = 'ocean_avg_';

file = struct('name',{});
file_len = 80;
%file_len = 40;



%fid = fopen(file_list);
%tline = fgetl(fid);
%while ischar(tline)
%    %disp(tline)
%    file_len = file_len + 1;
%    file(file_len).name = tline;
%    tline = fgetl(fid);
%end
%fclose(fid);

% Get grid dimensions
xi = length(ncread(gname,'xi_rho'));
eta = length(ncread(gname,'eta_rho'));

% Get grid
lat = ncread(gname,'lat_rho');
lon = ncread(gname,'lon_rho');
theta = ncread(gname,'angle');
mask_rho = ncread(gname,'mask_rho');
mask = mask_rho;
mask(mask == 0) = NaN;

% Empty array
curr_time = zeros(1,file_len);
ubar = zeros(xi,eta,file_len);
vbar = zeros(xi,eta,file_len);
sbar = zeros(xi,eta,file_len);
tbar = zeros(xi,eta,file_len);
NO3bar = zeros(xi,eta,file_len);
NH4bar = zeros(xi,eta,file_len);
phytoplanktonbar = zeros(xi,eta,file_len);
phytoplanktonNbar = zeros(xi,eta,file_len);

% Map projection
%m_proj('lambert','long',[8.0 13.2],'lat',[78.7 79.4]);
m_proj('lambert','long',[8.26 13.2],'lat',[78.76 79.38]);

% Define contour levels for temperature and salinity
contour_slines = 32:0.01:35.1;
%csrange = [ -Inf contour_slines ; contour_slines Inf ];
%csmap = parula(size(csrange,1));
%
contour_tlines = 3.0:0.05:5.0;
contour_NO3lines = 0.0:0.05:10.0;
contour_phytolines = 0:0.1:2.0;
contour_phytolinesN = 0:0.05:0.5;
%ctrange = [ -Inf contour_tlines ; contour_tlines Inf ];
%ctmap = summer(size(ctrange,1));

% Quiver density
q1 = 20:10:380;
q2 = 20:10:580;
qlen = length(q1)*length(q2);
qlon = [reshape(lon(q1,q2),[1 qlen]) 12.8];
qlat = [reshape(lat(q1,q2),[1 qlen]) 79.32];

% Set figure frame
Fig1 = figure('Position', [200 200 900 600]);
index = 0;
%index = 80;
for i=1:file_len
    index = index + 1; 
    %if i ~= 60
        
    if index < 10
      rfil = strcat(prename,'000',int2str(index),'.nc');
    end   
    if index >= 10 & i < 100
      rfil = strcat(prename,'00',int2str(index),'.nc');   
    end
    if index >=100
      rfil = strcat(prename,'0',int2str(index),'.nc');
    end  
    %
    
    file(i).name = rfil; 
    disp(file(i).name);
    fname = [dir_data file(i).name];
    %
    curr_time(i) = ncread(fname,'ocean_time')/86400.0 + t0;
    %
    if i==1
        % Find vertical depths
        fname
        Vtransform = ncread(fname,'Vtransform');
        Vstretching = ncread(fname,'Vstretching');
        theta_s = ncread(fname,'theta_s');
        theta_b = ncread(fname,'theta_b');
        hc = ncread(fname,'Tcline');
        h0 = ncread(fname,'h');
    end
    zeta = ncread(fname,'zeta');
    [z_w] = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 5, h0, zeta); %z_w = squeeze(abs(z_w));
    %z_w = depths(fname, gname, 5, 0, 1);   % ROMS function
    size(z_w)
    size(mask)
    for ii = 1:36
       z_w(:,:,ii) = z_w(:,:,ii) .* mask(:,:);     % mask land with NaN
    end
    %
    u = ncread(fname,'u');
    v = ncread(fname,'v');
    salt = ncread(fname,'salt');
    temp = ncread(fname,'temp');
    NO3 = ncread(fname,'NO3');
    NH4 = ncread(fname,'NH4');
    phytoplankton = ncread(fname,'phytoplankton');
    phytoplanktonN = ncread(fname,'Phy_N');
    %
    % Find velocity within depth range
    %k = find( (z_w > depth_range(1)) & (z_w < depth_range(2)) );
    for a=1:xi
        for b=1:eta
            if ~isnan(mask(a,b)) && z_w(a,b,1) < depth_range(2)
                % Depth of layer
                h_layer = 0;
                % Find first vertical index
                if depth_range(1) < z_w(a,b,1)
                    c = 1;
                else
                    c = find(z_w(a,b,:) < depth_range(1),1,'last');
                end
                % Depth integrate
                while (c < 36) && (z_w(a,b,c) < depth_range(2))
                    z0 = max( z_w(a,b,c), depth_range(1) );
                    z1 = min( z_w(a,b,c+1), depth_range(2) );
                    % Increase depth layer
                    h_layer = h_layer + (z1 - z0);
                    % Rotate velocity vector by angle
                    if (a==1) || (b==1) || (a==xi) || (b==eta)
                        ubar(a,b,i) = NaN;
                        vbar(a,b,i) = NaN;
                    else
                        U = ( u(a,b-1,c) + u(a,b,c) ) / 2.0;
                        V = ( v(a-1,b,c) + v(a,b,c) ) / 2.0;
                        ubar(a,b,i) = ubar(a,b,i) +  ...
                            ( U*cos( theta(a,b) ) -  ...
                              V*sin( theta(a,b) ) ) * (z1 - z0);
                        vbar(a,b,i) = vbar(a,b,i) +  ...
                            ( U*sin( theta(a,b) ) +  ...
                              V*cos( theta(a,b) ) ) * (z1 - z0);
                    end
                    % Salinity and Temperature
                    sbar(a,b,i) = sbar(a,b,i) + salt(a,b,c) * (z1 - z0);
                    tbar(a,b,i) = tbar(a,b,i) + temp(a,b,c) * (z1 - z0);
                    NO3bar(a,b,i) = NO3bar(a,b,i) + NO3(a,b,c) * (z1 - z0);
                    NH4bar(a,b,i) = NH4bar(a,b,i) + NH4(a,b,c) * (z1 - z0);
                    phytoplanktonbar(a,b,i) = phytoplanktonbar(a,b,i) + phytoplankton(a,b,c) * (z1 - z0);
                    phytoplanktonNbar(a,b,i) = phytoplanktonNbar(a,b,i) + phytoplanktonN(a,b,c) * (z1 - z0);
                    % increment c
                    c = c+1;
                end
                % Average over depth
                ubar(a,b,i) = ubar(a,b,i) / h_layer;
                vbar(a,b,i) = vbar(a,b,i) / h_layer;
                sbar(a,b,i) = sbar(a,b,i) / h_layer;
                tbar(a,b,i) = tbar(a,b,i) / h_layer;  
                NO3bar(a,b,i) = NO3bar(a,b,i) / h_layer;
                NH4bar(a,b,i) = NH4bar(a,b,i) / h_layer;
                phytoplanktonbar(a,b,i) = phytoplanktonbar(a,b,i) / h_layer;
                phytoplanktonNbar(a,b,i) = phytoplanktonNbar(a,b,i) / h_layer;
            else
                ubar(a,b,i) = NaN;
                vbar(a,b,i) = NaN;
                sbar(a,b,i) = NaN;
                tbar(a,b,i) = NaN;
                NO3bar(a,b,i) = NaN;
                NH4bar(a,b,i) = NaN;
                phytoplanktonbar(a,b,i) = NaN;
                phytoplanktonNbar(a,b,i) = NaN;
            end
        end
    end
    %
    % Plot daily average
    %
    m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
    hold on;
    %-- Salinity
%     [CS,CH] = m_contourf(lon,lat,sbar(:,:,i), ...
%         contour_slines,'LineStyle','none');
%     CMH = colormap(cm_haline);
    %-- Temperature
    [CS,CH] = m_contourf(lon,lat,tbar(:,:,i), ...
        contour_tlines,'LineStyle','none');
    CMH = colormap(cm_thermal);
    %--
    qubar = [reshape(ubar(q1,q2,i),[1 qlen]) 0.1];
    qvbar = [reshape(vbar(q1,q2,i),[1 qlen]) 0.0];
    m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
    m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
    %--
    hold off;
    m_grid('box','fancy','tickdir','in');
    set(gca,'Color',CMH(1,:));
    m_contfbar(1.03, [0.1 0.9], CS, CH, 'endpiece','yes');
    pause(0.1);
    m_contfbar();
    %end
end

%save(sprintf('vel_salt_temp_%s.mat', label), 'ubar', 'vbar', 'sbar', 'tbar', 'curr_time');

%
% Monthly mean values
%
ubar_mean = mean(ubar,3);
vbar_mean = mean(vbar,3);
sbar_mean = mean(sbar,3);
tbar_mean = mean(tbar,3);
NO3_mean = mean(NO3bar,3);
NH4_mean = mean(NH4bar,3);
phytoplankton_mean = mean(phytoplanktonbar,3);
phytoplanktonN_mean = mean(phytoplanktonNbar,3);

% Quiver
qubar = [reshape(ubar_mean(q1,q2),[1 qlen]) 0.1];
qvbar = [reshape(vbar_mean(q1,q2),[1 qlen]) 0.0];

%
% Mean salinity
%
Fig2 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- Salinity
[CS,CH] = m_contourf(lon,lat,sbar_mean,  ...
    contour_slines,'LineStyle','none');
CMH = colormap(cm_haline);
%-- Quiver
m_quiver(qlon,qlat,qubar,qvbar,'Color','red');
m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMH(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH, 'endpiece','yes');
title(sprintf('Salinity, %s %s %s', case_label,month,range_label));
%print('-dpng','-r300',sprintf('salt_%s.png', label));

%
% Mean temperature
%
Fig3 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- Temperature
[CS,CH] = m_contourf(lon,lat,tbar_mean, ...
    contour_tlines,'LineStyle','none');
CMT = colormap(cm_thermal);
%-- Quiver
m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
%m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','yes');
title(sprintf('Temperature, %s %s %s', case_label,month,range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig4 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- NO3
[CS,CH] = m_contourf(lon,lat,NO3_mean, ...
    contour_NO3lines,'LineStyle','none');
CMT = colormap(cm_thermal);
%-- Quiver
m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
%m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','yes');
title(sprintf('NO3, %s %s %s', case_label,month,range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig5 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- NH4
[CS,CH] = m_contourf(lon,lat,NH4_mean, ...
    contour_NO3lines,'LineStyle','none');
CMT = colormap(cm_thermal);
%-- Quiver
m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
%m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','yes');
title(sprintf('NH4, %s %s %s', case_label,month,range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig6 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- phytoplankton
[CS,CH] = m_contourf(lon,lat,phytoplankton_mean, ...
    contour_phytolines,'LineStyle','none');
CMT = colormap(cm_thermal);
%-- Quiver
m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
%m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','yes');
title(sprintf('Phytoplankton, %s %s %s', case_label,month,range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig7 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- phytoplanktonN
[CS,CH] = m_contourf(lon,lat,phytoplanktonN_mean, ...
    contour_phytolinesN,'LineStyle','none');
CMT = colormap(cm_thermal);
%-- Quiver
m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
%m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','yes');
title(sprintf('PhytoplanktonN, %s %s %s', case_label,month,range_label));