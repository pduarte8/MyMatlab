function [NPP_mean,lon,lat,mask_rho] = Zlevel_net_production(NPParray)


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

%diag_file = strcat(dir_data,'ocean_dia.nc')
%NPP = ncread(diag_file,'P_Production');
NPP = NPParray;
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

NPPbar = zeros(xi,eta,file_len);

% Map projection
%m_proj('lambert','long',[8.0 13.2],'lat',[78.7 79.4]);
m_proj('lambert','long',[8.26 13.2],'lat',[78.76 79.38]);

% Define contour levels for temperature and salinity


contour_NPPlines = -0.00005:0.000005:0.0005;


% Quiver density
q1 = 20:10:380;
q2 = 20:10:580;
qlen = length(q1)*length(q2);
qlon = [reshape(lon(q1,q2),[1 qlen]) 12.8];
qlat = [reshape(lat(q1,q2),[1 qlen]) 79.32];

% Set figure frame

index = 0;
%index = 80;

for i=1:file_len
    index = index + 1; 
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
                    
                    % Production
                    
                    NPPbar(a,b,i) = NPPbar(a,b,i) + NPP(a,b,c,index) * (z1 - z0);
                    % increment c
                    c = c+1;
                end
                % Average over depth               
                NPPbar(a,b,i) = NPPbar(a,b,i) / h_layer;
            else
                NPPbar(a,b,i) = NaN;
            end
        end
    end
end

%save(sprintf('vel_salt_temp_%s.mat', label), 'ubar', 'vbar', 'sbar', 'tbar', 'curr_time');

%
% Monthly mean values

NPP_mean = mean(NPPbar,3);



Fig1 = figure('Position', [200 200 900 600]);
hold on;
%-- Land contour
%m_contour(lon,lat,mask_rho,'LineWidth',2.0,'LineColor','black');
%-- NPP
[CS,CH] = m_contourf(lon,lat,NPP_mean, ...
    contour_NPPlines,'LineStyle','none');
CMT = colormap(cm_thermal);
hold off;
m_grid('box','fancy','tickdir','in');
set(gca,'Color',CMT(1,:));
%m_contfbar(1.03, [0.1 0.9], CS, CH,  'endpiece','yes');
title(sprintf('NPP, %s %s %s', case_label,month,range_label));
%print('-dpng','-r300',sprintf('temp_%s.png', label));

