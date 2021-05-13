clear;

% Plot parameters
sim_label = 'K160v2';
%
% case_label = 'subglacial';
% case_label = 'surface';
% case_label = 'future';
%
 month = 'June';
% month = 'July';
% month = 'Aug';
%month = 'Sep';

% Depth range : water depths are negative
%               use large positive value e.g. 99 to include sea surface
 depth_range = [ -30   99];  range_label = '0 to 30m';
% depth_range = [ -30  -20];  range_label = '20 to 30m';
% depth_range = [ -55  -45];  range_label = '45 to 55m';
% depth_range = [ -80  -70];  range_label = '70 to 80m';
% depth_range = [-105  -95];  range_label = '95 to 105m';
% depth_range = [ -10   1000];  range_label = '0 to 10m'; %The one I used
% depth_range = [ -5   1000];  range_label = '0 to 5m'; 
% depth_range = [ -30   1000];  range_label = '0 to 30m'; 



t0 = datenum([1948 1 1 0 0 0]);

N = 35;
for expno=1:3
switch expno
case 1; case_label = 'surface';
case 2; case_label = 'subglacial';
case 3; case_label ='future';    
end    
% Grid file
switch case_label
    case 'subglacial'
        % Grid file
        gname = '/cluster/work/users/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/kongsfjorden_160m_grid_present.nc';
        % Data file directory
        dir_data = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/Sim14/';
        % file list
        file_list = sprintf('files_%s_%s_%s.txt', sim_label,case_label,month);
    case 'future'
        % Grid file
        gname = '/cluster/work/users/pduarte/tmproms/run/run_Kongsfjorden-160m_future/kongsfjorden_160m_grid_future.nc';
        % Data file directory
        dir_data = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/Sim18/';
        %dir_data = '/cluster/work/users/pduarte/tmproms/run/run_Kongsfjorden-160m_future/';
        % file list
        file_list = sprintf('files_%s_%s_%s.txt', sim_label,case_label,month);
    case 'surface'
        % Grid file
        gname = '/cluster/work/users/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/kongsfjorden_160m_grid_present.nc';
        % Data file directory
        dir_data = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/Sim16/';
        % file list
        file_list = sprintf('files_%s_subglacial_%s.txt', sim_label,month);
    otherwise
        error('Not a valid case.')
end
label = sprintf('%s_%s_%s', case_label,month,range_label);
prename = 'ocean_avg_';

file = struct('name',{});

file_len = 80;
%file_len = 40;

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
phytoplanktonbar = zeros(xi,eta,file_len);

size(sbar)
% Map projection
m_proj('lambert','long',[11.0 13.2],'lat',[78.76 79.2]);
%m_proj('lambert','long',[8.26 13.2],'lat',[78.76 79.38]);

% Define contour levels for velocity
contour_vlines = -0.2:0.1:0.2;
% Define contour levels for temperature and salinity
contour_slines = -1:0.1:1.0;
%csrange = [ -Inf contour_slines ; contour_slines Inf ];
%csmap = parula(size(csrange,1));
%
contour_tlines = -1:0.1:1.0;
contour_NO3lines = -2:0.1:2.0;
contour_phytolines = -1.0:0.1:1.0;
%ctrange = [ -Inf contour_tlines ; contour_tlines Inf ];
%ctmap = summer(size(ctrange,1));

% Quiver density
q1 = 1:4:402;
q2 = 1:4:602;
qlen = length(q1)*length(q2);
qlon = [reshape(lon(q1,q2),[1 qlen]) 12.8];
qlat = [reshape(lat(q1,q2),[1 qlen]) 79.32];
qlon(1) = 12.8;
qlat(1) = 79.05;


index = 0;
%index = 80;
for i=1:file_len
    index = index + 1; 
    %if index ~= 60
        
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
    phytoplankton = ncread(fname,'phytoplankton');
    
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
                    phytoplanktonbar(a,b,i) = phytoplanktonbar(a,b,i) + phytoplankton(a,b,c) * (z1 - z0);
                    % increment c
                    c = c+1;
                end
                % Average over depth
                ubar(a,b,i) = ubar(a,b,i) / h_layer;
                vbar(a,b,i) = vbar(a,b,i) / h_layer;
                sbar(a,b,i) = sbar(a,b,i) / h_layer;
                tbar(a,b,i) = tbar(a,b,i) / h_layer;  
                NO3bar(a,b,i) = NO3bar(a,b,i) / h_layer;
                phytoplanktonbar(a,b,i) = phytoplanktonbar(a,b,i) / h_layer;
            else
                ubar(a,b,i) = NaN;
                vbar(a,b,i) = NaN;
                sbar(a,b,i) = NaN;
                tbar(a,b,i) = NaN;
                NO3bar(a,b,i) = NaN;
                phytoplanktonbar(a,b,i) = NaN;
            end
        end
    end
    qubar = [reshape(ubar(q1,q2,i),[1 qlen]) 0.1];
    qvbar = [reshape(vbar(q1,q2,i),[1 qlen]) 0.0];
end



%
% Monthly mean values
%
ubar_mean = mean(ubar,3);
vbar_mean = mean(vbar,3);
for ii = 1:402
     for jj = 1:602
        res_vel_bar(ii,jj) = sqrt(ubar_mean(ii,jj)^2+vbar_mean(ii,jj)^2) * mask(ii,jj);
     end
end     
sbar_mean = mean(sbar,3);
tbar_mean = mean(tbar,3);
NO3_mean = mean(NO3bar,3);
phytoplankton_mean = mean(phytoplanktonbar,3);

% Quiver
qubar = [reshape(ubar_mean(q1,q2),[1 qlen]) 0.1];
qvbar = [reshape(vbar_mean(q1,q2),[1 qlen]) 0.0];
switch case_label
    case 'surface'; 
        mask_surf = mask_rho;
        res_vel_bar_surf = res_vel_bar;
        sbar_mean_surf = sbar_mean;
        tbar_mean_surf = tbar_mean;
        NO3_mean_surf = NO3_mean;
        phytoplankton_mean_surf = phytoplankton_mean;
        qubar_surf = qubar;
        qvbar_surf = qvbar;
    case 'subglacial'; 
        mask_sub = mask_rho;
        res_vel_bar_sub = res_vel_bar;
        sbar_mean_sub = sbar_mean;
        tbar_mean_sub = tbar_mean;
        NO3_mean_sub = NO3_mean;
        phytoplankton_mean_sub = phytoplankton_mean;
        qubar_sub = qubar;
        qvbar_sub = qvbar;
    case 'future'; 
        mask_fut = mask_rho;
        res_vel_bar_fut = res_vel_bar;
        sbar_mean_fut = sbar_mean;
        tbar_mean_fut = tbar_mean;
        NO3_mean_fut = NO3_mean;
        phytoplankton_mean_fut = phytoplankton_mean;
        qubar_fut = qubar;
        qvbar_fut = qvbar;  
end 
case_label
end %cases
%Mean velocity

Fig1 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_sub,'LineWidth',2.0,'LineColor','black');
%-- Velocity
diff = res_vel_bar_surf-res_vel_bar_sub;
diff = diff.*mask_sub;
[CS,CH] = m_contourf(lon,lat,diff,  ...
    contour_vlines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMH = colormap(cm_haline);
%-- Quiver

%m_quiver(qlon,qlat,qubar,qvbar,'Color','red');
qubar_sub(1) = 0.1;
qvbar_sub(1) = 0.0;
m_quiver(qlon(1:6350),qlat(1:6350),qubar_sub(1:6350),qvbar_sub(1:6350),1.5,'Color','red');
m_text(qlon(1)-0.05,qlat(1)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMH(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH, 'endpiece','no');
title(sprintf('Current speed, %s %s %s', 'PSD-PTG',range_label));

Fig2 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_fut,'LineWidth',2.0,'LineColor','black');
%-- Velocity
diff = res_vel_bar_fut-res_vel_bar_sub;
diff = diff.*mask_sub;
[CS,CH] = m_contourf(lon,lat,diff,  ...
    contour_vlines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_fut-mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMH = colormap(cm_haline);
%-- Quiver

%m_quiver(qlon,qlat,qubar,qvbar,'Color','red');
qubar_sub(1) = 0.1;
qvbar_sub(1) = 0.0;
m_quiver(qlon(1:6350),qlat(1:6350),qubar_sub(1:6350),qvbar_sub(1:6350),1.5,'Color','red');
m_text(qlon(1)-0.05,qlat(1)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMH(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH, 'endpiece','no');
title(sprintf('Current speed, %s %s %s', 'FLG-PTG',range_label));
% Mean salinity
%
Fig3 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_sub,'LineWidth',2.0,'LineColor','black');
m_contour(lon,lat,h0.*(mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
%-- Salinity
[CS,CH] = m_contourf(lon,lat,(sbar_mean_surf-sbar_mean_sub).*mask_sub,  ...
    contour_slines,'LineStyle','none');
CMH = colormap(cm_haline);
%-- Quiver
%m_quiver(qlon,qlat,qubar,qvbar,'Color','red');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMH(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH, 'endpiece','no');
title(sprintf('Salinity, %s %s %s', 'PSD-PTG',range_label));
%print('-dpng','-r300',sprintf('salt_%s.png', label));

Fig4 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
m_contour(lon,lat,h0.*(mask_fut-mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
%-- Salinity
[CS,CH] = m_contourf(lon,lat,(sbar_mean_fut-sbar_mean_sub).*mask_sub,  ...
    contour_slines,'LineStyle','none');
CMH = colormap(cm_haline);
%-- Quiver
%m_quiver(qlon,qlat,qubar,qvbar,'Color','red');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMH(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH, 'endpiece','no');
title(sprintf('Salinity, %s %s %s', 'FLG-PTG',range_label));
%print('-dpng','-r300',sprintf('salt_%s.png', label));
%
% Mean temperature
%
Fig5 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_sub,'LineWidth',2.0,'LineColor','black');
%-- Temperature
[CS,CH] = m_contourf(lon,lat,(tbar_mean_surf-tbar_mean_sub).*mask_sub, ...
    contour_tlines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMT = colormap(cm_thermal);
%-- Quiver
%m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','no');
title(sprintf('Temperature, %s %s %s', 'PSD-PTG',range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig6 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- Temperature
[CS,CH] = m_contourf(lon,lat,(tbar_mean_fut-tbar_mean_sub).*mask_sub, ...
    contour_tlines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_fut-mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMT = colormap(cm_thermal);
%-- Quiver
%m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','no');
title(sprintf('Temperature, %s %s %s', 'FLG-PTG',range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig7 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_sub,'LineWidth',2.0,'LineColor','black');
%-- NO3
[CS,CH] = m_contourf(lon,lat,(NO3_mean_surf-NO3_mean_sub).*mask_sub, ...
    contour_NO3lines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMT = colormap(cm_thermal);
%-- Quiver
%m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','no');
title(sprintf('NO3, %s %s %s','PSD-PTG',range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig8 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- NO3
[CS,CH] = m_contourf(lon,lat,(NO3_mean_fut-NO3_mean_sub).*mask_sub, ...
    contour_NO3lines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_fut-mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMT = colormap(cm_thermal);
%-- Quiver
%m_quiver(qlon,qlat,qubar,qvbar,'Color','cyan');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','no');
title(sprintf('NO3, %s %s %s','FLG-PTG',range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig9 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_sub,'LineWidth',2.0,'LineColor','black');
%-- phytoplankton
[CS,CH] = m_contourf(lon,lat,(phytoplankton_mean_surf-phytoplankton_mean_sub).*mask_sub, ...
    contour_phytolines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMT = colormap(cm_thermal);
%-- Quiver
%m_quiver(qlon(1:3900),qlat(1:3900),qubar(1:3900),qvbar(1:3900),'Color','cyan');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','no');
title(sprintf('Phytoplankton, %s %s %s', 'PSD-PTG',month,range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

Fig10 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- phytoplankton
[CS,CH] = m_contourf(lon,lat,(phytoplankton_mean_fut-phytoplankton_mean_sub).*mask_sub, ...
    contour_phytolines,'LineStyle','none');
m_contour(lon,lat,h0.*(mask_fut-mask_sub),[0 5 10 15 25 50 75 100 125 150 175 200 250 300 400 500 1000],'LineWidth',1.0,'LineColor','black');
CMT = colormap(cm_thermal);
%-- Quiver
%m_quiver(qlon(1:3900),qlat(1:3900),qubar(1:3900),qvbar(1:3900),'Color','cyan');
%m_text(qlon(end)-0.05,qlat(end)+0.015,'10 cm/s');
%--
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','no');
title(sprintf('Phytoplankton, %s %s %s', 'FLG-PTG',range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));
