clear;

% Plot parameters
sim_label = 'K160_Sim25';
%
case_label = 'subglacial';
% case_label = 'surface';
%case_label = 'future';
%
 month = 'June';
% month = 'July';
% month = 'Aug';
%month = 'Sep';

% Depth range : water depths are negative
%               use large positive value e.g. 99 to include sea surface
%depth_range = [ -10   99];  range_label = '0to10m';
 depth_range = [ -30  -20];  range_label = '20to30m';
% depth_range = [ -55  -45];  range_label = '45to55m';
% depth_range = [ -80  -70];  range_label = '70to80m';
% depth_range = [-105  -95];  range_label = '95to105m';

label = sprintf('%s_%s_%s', case_label,month,range_label);

t0 = datenum([1948 1 1 0 0 0]);

% Grid file
switch case_label
    case 'subglacial'
        % Grid file
        gname = '/cluster/shared/arcticfjord/run_Kongsfjorden-160m_present_subglacial/kongsfjorden_160m_grid_present.nc';
        % Data file directory
        dir_data = '/cluster/work/users/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/';
        %dir_data = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/SecondRoundOfSimulations/Sim23_repeated/';
        % file list
        %file_list = sprintf('files_%s_%s_%s.txt', sim_label,case_label,month);
    case 'future'
        % Grid file
        gname = 'D:\Data\ROMS\norkyst\NorFjords-160m_Forcing\Kongsfjorden-160m_v2_future\Grid\kongsfjorden_160m_grid_future.nc';
        % Data file directory
        dir_data = 'D:\Data\roms_run\Kongsfjorden-160m\K160_v2\future\daily_avg\';
        % file list
        %file_list = sprintf('files_%s_%s_%s.txt', sim_label,case_label,month);
    case 'surface'
        % Grid file
        gname = 'D:\Data\ROMS\norkyst\NorFjords-160m_Forcing\Kongsfjorden-160m_v2_present\Grid\kongsfjorden_160m_grid_present_NEW.nc';
        % Data file directory
        dir_data = 'D:\Data\roms_run\Kongsfjorden-160m\K160_v2\present\daily_avg\';
        % file list
        %file_list = sprintf('files_%s_subglacial_%s.txt', sim_label,month);
    otherwise
        error('Not a valid case.')
end

prename = 'ocean_avg';

file = struct('name',{});
file_len = 1;

%fid = fopen(file_list);
%tline = fgetl(fid);
%while ischar(tline)
    %disp(tline)
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

% Map projection
%m_proj('lambert','long',[8.0 13.2],'lat',[78.7 79.4]);
m_proj('lambert','long',[8.26 13.2],'lat',[78.76 79.38]);

% Define contour levels for temperature and salinity
contour_slines = 31.0:0.1:35.0;
%csrange = [ -Inf contour_slines ; contour_slines Inf ];
%csmap = parula(size(csrange,1));
%
contour_tlines = 2.0:0.1:6.0;
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
index = 36;
for i=1:file_len
    %disp(file(i).name);
    index = index + 1;
    if index < 10
      rfil = strcat(prename,'000',int2str(index),'.nc');
    end   
    if index >= 10 & index < 100
      rfil = strcat(prename,'00',int2str(index),'.nc');   
    end
    if index >=100
      rfil = strcat(prename,'0',int2str(index),'.nc');
    end  
    %
    
    file(i).name = rfil;
    
    fname = [dir_data file(i).name]
    
    curr_time(i) = ncread(fname,'ocean_time')/86400.0 + t0
    %
    z_w = depths(fname, gname, 5, 0, 1);   % ROMS function
    z_w = z_w(:,:,1:end) .* mask(:,:);     % mask land with NaN
    %
    u = ncread(fname,'u');
    v = ncread(fname,'v');
    salt = ncread(fname,'salt');
    temp = ncread(fname,'temp');
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
                    if (a==1) || (b==1) || (a==xi) || (b==eta) ||...
                       length(find(u(a,b,:) > 10) > 0) || length(find(v(a,b,:) > 10) > 0) ||...
                       ((b > 0) && length(find(u(a,b-1,:) > 10) > 0)) ||...
                       ((a > 0) && length(find(v(a-1,b,:) > 10) > 0))
                        ubar(a,b,i) = NaN;
                        vbar(a,b,i) = NaN;
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
                    % increment c
                    c = c+1;
                end
                % Average over depth
                ubar(a,b,i) = ubar(a,b,i) / h_layer;
                vbar(a,b,i) = vbar(a,b,i) / h_layer;
                sbar(a,b,i) = sbar(a,b,i) / h_layer;
                tbar(a,b,i) = tbar(a,b,i) / h_layer;                
            else
                ubar(a,b,i) = NaN;
                vbar(a,b,i) = NaN;
                sbar(a,b,i) = NaN;
                tbar(a,b,i) = NaN;
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
end

save(sprintf('vel_salt_temp_%s.mat', label), 'ubar', 'vbar', 'sbar', 'tbar', 'curr_time');

%
% Monthly mean values
%
ubar_mean = mean(ubar,3);
vbar_mean = mean(vbar,3);
sbar_mean = mean(sbar,3);
tbar_mean = mean(tbar,3);

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
print('-dpng','-r300',sprintf('salt_%s.png', label));

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
m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','yes');
title(sprintf('Temperature, %s %s %s', case_label,month,range_label));
print('-dpng','-r300',sprintf('temp_%s.png', label));
