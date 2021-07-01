clear all; close all;
%vid = VideoWriter('velocity_animation.avi','Motion JPEG AVI');
%vid.Quality = 100;
%vid.FrameRate = 15;
%open(vid);
%inPath = check_path('/home/pduarte/matlabmfiles/ROMS');
%if (~inPath)
%    disp('Add ROMS toolbox to path');
    run('/cluster/home/pduarte/matlabmfiles/ROMS/matlab/startup.m');
%end

% Calculate and plot mean current speed through a defined section
% Calculate volume flux through a defined section according to depths

% Transect denoted with grid points (must be parallel to grid axis) (EAST)
% i0 =  82; i1 = 105; j0 =  81; j1 =  73; MapFile = 'Kronebreen_transect';
% i0 = 113; i1 = 134; j0 = 144; j1 = 132; MapFile ='Ny-Ålesund transect';
% i0 = 143; i1 = 156; j0 =  76; j1 =  76; MapFile = 'Kongsbreen_transect';
% i0 = 138; i1 = 201; j0 = 214; j1 = 180; MapFile = 'Kapp_Guissez_transect'; 
% i0 =  90; i1 = 107; j0 =  91; j1 =  91; MapFile = 'Lovenøyane_transect';
%i0 =  2; i1 = 400; j0 =  600; j1 = 600; MapFile = 'W_boundary_transect';
%i0 =  2; i1 = 2; j0 =  318; j1 = 600; MapFile = 'S_boundary_transect';
%i0 =  2; i1 = 2; j0 =  176; j1 = 244; MapFile = 'S_inner_boundary_transect';
% i0 = 132; i1 = 177; j0 = 174; j1 = 150; MapFile = 'Kongsfjorden_transect_K160v2.png';
i0 =  400; i1 = 400; j0 =  278; j1 = 600; MapFile = 'N_boundary_transect';
%rdir = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_present_subglacial/Sim1'
%rdir = '/global/work/pduarte/tmproms/run/run_Kongsfjorden-160m_future/Sim2/StalloAndFramResults'
%rdir = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/Sim1'
%rdir = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/Sim5'
rdir = '/tos-project3/NS9081K/NORSTORE_OSL_DISK/NS9081K/K160_bgc/ThirdRoundOfSimulations/Sim31'
prename = 'ocean_avg_';
nfiles = 5;
%nfiles = 40;
t0 = datenum([2007 06 01 0 0 0]);  % Last day before beginning of data reco
%t0 = datenum([2007 08 20 0 0 0]);  % Last day before beginning of data record
z1 = -10; z2 = 500;  % Upper and lower depth (positive numbers in meter, z1 < z2). If all depths: make sure z2 is large enough

% Start loop over files
first = 1;  % Logical switch to prevent unnecessary repetition of actions
%for fn=1:nfiles
index = 0;
%index = 80; 
aindex = index;
for fn = 1:nfiles
  index = index + 1;  
  if index < 10
     rfil = strcat(rdir,'/',prename,'000',int2str(index),'.nc');
  end   
  if index >= 10 & fn < 100
     rfil = strcat(rdir,'/',prename,'00',int2str(index),'.nc');   
  end
  if index >=100
     rfil = strcat(rdir,'/',prename,'0',int2str(index),'.nc');
  end  
  
  if 0  % To find transect based on average (surface) currents
    mask = ncread(rfil,'mask_rho'); mask(mask==0) = NaN;
    h = ncread(rfil,'h'); h=h.*mask;
    figure(2); clf; contourf(h'); colorbar; hold on;
                    line([i0,i1],[j0,j1],'Color','r');
                    axis equal;
                    drawnow;
                    %print('-dpng','-r600',MapFile);
    clear mask h
  end  % if 0
  if first

    % Read constant model fields
    mask = ncread(rfil,'mask_rho');
    lon = ncread(rfil,'lon_rho');
    lat = ncread(rfil,'lat_rho');
    phi = ncread(rfil,'angle');
    %pm = ncread(rfil,'pm'); dx = 1./pm; clear pm

    % Define line X,Y from (i0,j0) to (i1,j1) with spacing according to model resolution
    k = 1; X(k) = i0; Y(k) = j0;
    DX = i1 - i0; if DX ~= 0; signx = DX/abs(DX); else signx = 1; end
    DY = j1 - j0; if DY ~= 0; signy = DY/abs(DY); else signy = 1; end
    if abs(DX) >= abs(DY)
      xtrans = true;
      for j=1:abs(DX)
        k = k + 1;
        X(k) = i0 + signx*j;
        Y(k) = j0 + signy*j*abs(DY)/abs(DX);
      end  % j
    else
      xtrans = false;
      for j=1:abs(DY)
        k = k + 1;
        X(k) = i0 + signx*j*abs(DX)/abs(DY);
        Y(k) = j0 + signy*j;
      end  % j
    end
    X = round(X); Y = round(Y);
    clear k DX DY j 
    
    % Find horizontal discretization along transect
    if xtrans  % use U-edge referenc points
        lon_u = ncread(rfil,'lon_u');
        lat_u = ncread(rfil,'lat_u');
        if signx > 0
            rng_lon(1) = lon_u(i0,j0);
            rng_lon(2) = lon_u(i1+1,j1);
            rng_lat(1) = lat_u(i0,j0);
            rng_lat(2) = lat_u(i1+1,j1);
        else
            rng_lon(1) = lon_u(i1,j1);
            rng_lon(2) = lon_u(i0+1,j0);
            rng_lat(1) = lat_u(i1,j1);
            rng_lat(2) = lat_u(i0+1,j0);
        end
        clear lon_u lat_u
    else
        lon_v = ncread(rfil,'lon_v');
        lat_v = ncread(rfil,'lat_v');
        if signy > 0
            rng_lon(1) = lon_v(i0,j0);
            rng_lon(2) = lon_v(i1,j1+1);
            rng_lat(1) = lat_v(i0,j0);
            rng_lat(2) = lat_v(i1,j1+1);
        else
            rng_lon(1) = lon_v(i1,j1);
            rng_lon(2) = lon_v(i0,j0+1);
            rng_lat(1) = lat_v(i1,j1);
            rng_lat(2) = lat_v(i0,j0+1);
        end
        clear lon_v lat_v
    end
    clear i0 i1 j0 j1 %signx signy

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

    % Pick out values of mask, lon, lat, phi, h and z along transect X,Y
    for i=1:length(X)
      MASK(i) = mask(X(i),Y(i));
      LON(i) = lon(X(i),Y(i));
      LAT(i) = lat(X(i),Y(i));
      PHI(i) = phi(X(i),Y(i));
      H(i) = h0(X(i),Y(i));
      Z0(i,:) = z0(X(i),Y(i),:);
      Z0_W(i,:) = z0_w(X(i),Y(i),:);
    end
    clear mask lon lat phi z0 z0_w
    mask = MASK; lon = LON; lat = LAT; phi = PHI; h = H; z0 = Z0; z0_w = Z0_W;
    clear MASK LON LAT PHI H Z0 Z0_W

    % Find actual distance between points
    rng = m_lldist(rng_lon,rng_lat)*1000.;
    for i=1:length(X)
        dx(i) = rng/length(X);
    end
    
%     % Find actual distance between points
%     range = m_lldist(lon,lat)*1000.;
%     for i=1:length(lon)
%       if i == 1
%         dx(i) = range(i);
%       elseif i == length(lon)
%         dx(i) = range(i-1);
%       else
%         dx(i) = 0.5*(range(i-1) + range(i));
%       end
%     end
%     clear range

    % Find angle based on the orientation of the transect
    a1 = max(lat) - min(lat);
    a2 = (max(lon) - min(lon))*cosd(mean(lat));
    angle_transect = atan2(a1,a2);
    clear a1 a2

    % Find vertical spacing
    for i=1:length(lon)
      for k=1:N
        dz0(i,k) = z0_w(i,k) - z0_w(i,k+1);
        %if k == 1
        %  dz(i,k) = h(i) - 0.5*(z(i,k) + z(i,k+1));
        %elseif k == N
        %  dz(i,k) = 0.5*(z(i,k) + z(i,k-1));
        %else
        %  dz(i,k) = 0.5*(z(i,k-1) - z(i,k+1));
        %end
      end  % k
      % Set dz0=0 if above or below upper and/or lower depth
      for k=1:N
        if z0(i,k) > z2 | z0(i,k) < z1
          dz0(i,k) = 0.;
        end
      end
      % Adjust dz0 according to z2 and z1
      k = find(dz0(i,:) > 0);
      if k(1) > 1;   dz0(i,k(1)) = z2 - 0.5*(z0(i,k(1)) + z0(i,k(1)+1)); end
      if k(end) < N; dz0(i,k(end)) = 0.5*(z0(i,k(end)) + z0(i,k(end)-1)) - z1; end
      clear k
    end  % i

    % Initialize u_transect
    %z_transect = zeros(length(X),N,nfiles);
    zeta_transect = zeros(length(X),nfiles);
    dz_transect = zeros(length(X),N,nfiles);
    %u_transect = zeros(length(X),N,nfiles);
    %temp_transect = zeros(length(X),N,nfiles);
    %salt_transect = zeros(length(X),N,nfiles);
    %no3_transect = zeros(length(X),N,nfiles);
    %phyto_transect = zeros(length(X),N,nfiles);

    first = 0;
  end  % if first

  tid = ncread(rfil,'ocean_time'); tid = tid/(24*3600); dato = gregorian(tid + julian(1948,1,1));

  % Find vertical depths
  %zeta = 0.;
  zeta = ncread(rfil,'zeta');
  [z] = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 1, h0, zeta); z = squeeze(abs(z));
  [z_w] = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, 5, h0, zeta); z_w = squeeze(abs(z_w));
  for i=1:length(X)
     Z(i,:) = z(X(i),Y(i),:);
     Z_W(i,:) = z_w(X(i),Y(i),:);
  end
  zeta_transect(:,fn) = squeeze(Z_W(:,end));
  clear z z_w

  % Find vertical spacing
  for i=1:length(lon)
    for k=1:N
      if Z(i,k) > z2 | Z(i,k) < z1
        dz_transect(i,k,fn) = 0.;
      else
        dz_transect(i,k,fn) = Z_W(i,k) - Z_W(i,k+1);
      end
    end  % k
    % Adjust dz0 according to z2 and z1
    k = find(dz_transect(i,:,fn) > 0);
    if k(1) > 1;   dz_transect(i,k(1),fn) = z2 - 0.5*(Z(i,k(1)) + Z(i,k(1)+1)); end
    if k(end) < N; dz_transect(i,k(end),fn) = 0.5*(Z(i,k(end)) + Z(i,k(end)-1)) - z1; end
    clear k
  end
  clear Z Z_W
   
  % Read current vectors and extract velocities along transect
  uin = ncread(rfil,'u');
  vin = ncread(rfil,'v');
  temp = ncread(rfil,'temp');
  salt = ncread(rfil,'salt');
  %no3 = ncread(rfil,'NO3');
  %phyto = ncread(rfil,'phytoplankton');
  
  
  for i=1:length(X)
    for j=1:N
%      UIN(i,j) = nanmean([uin(X(i),Y(i),j) uin(X(i)+1,Y(i),j)]);
%      VIN(i,j) = nanmean([vin(X(i),Y(i),j) vin(X(i),Y(i)+1,j)]);
      UIN(i,j) = nanmean([uin(X(i)-1,Y(i),j) uin(X(i),Y(i),j)]);
      VIN(i,j) = nanmean([vin(X(i),Y(i)-1,j) vin(X(i),Y(i),j)]);
    end
%     UIN(i,:) = uin(X(i),Y(i),:);
%     VIN(i,:) = vin(X(i),Y(i),:);
    TEMP(i,:) = temp(X(i),Y(i),:);
    SALT(i,:) = salt(X(i),Y(i),:);
    %NO3(i,:) = no3(X(i),Y(i),:);
    %PHYTO(i,:) = phyto(X(i),Y(i),:);
  end
  UIN(isnan(UIN)) = 0.;           % Avoid setting grid cell velocity to
  VIN(isnan(VIN)) = 0.;           % NaN if one component is finite
  clear uin vin
  uin = UIN; vin = VIN;
  temp = TEMP; salt = SALT; %no3 = NO3; phyto = PHYTO; 
  clear UIN VIN TEMP SALT NO3 PHYTO 

  % Rotate vectors according to east-west and north-south
  for i=1:length(X)
    for k=1:N
      uin_EW(i,k) = uin(i,k)*cos(phi(i)) - vin(i,k)*sin(phi(i));
      vin_NS(i,k) = vin(i,k)*cos(phi(i)) + uin(i,k)*sin(phi(i));
    end
  end

  % Rotate vectors according to defined transect
  % NB: v_normal is now directed normal to transect AND positive values mean current OUT OF the fjord
  for i=1:length(X)
    for k=1:N
      u_along(i,k) = uin_EW(i,k)*cos(angle_transect) + vin_NS(i,k)*sin(angle_transect);
      v_normal(i,k) = vin_NS(i,k)*cos(angle_transect) - uin_EW(i,k)*sin(angle_transect);
    end
  end

  if 0  % Plot to test that rotation of vectors is ok
    lon_rho = ncread(rfil,'lon_rho'); lon_rho = lon_rho(2:end-1,2:end-1);
    lat_rho = ncread(rfil,'lat_rho'); lat_rho = lat_rho(2:end-1,2:end-1);
    mask_rho = ncread(rfil,'mask_rho'); mask_rho = mask_rho(2:end-1,2:end-1);
    u = ncread(rfil,'u'); u = squeeze(u(:,:,end)); u = 0.5*(u(1:end-1,2:end-1) + u(2:end,2:end-1));
    v = ncread(rfil,'v'); v = squeeze(v(:,:,end)); v = 0.5*(v(2:end-1,1:end-1) + v(2:end-1,2:end));
    figure(1); clf; quiver(u',v',5.,'k-','LineWidth',0.2); hold on;
                    contour(mask_rho',[0.5 0.5],'k'); hold on;
                    plot(X,Y,'r*');
                    axis equal; axis([360 375 25 50]);
    figure(2); clf; m_proj('lambert','long',[12.35 12.50],'lat',[78.88 78.92]);
                    m_quiver(lon,lat,uin_EW(:,end)',vin_NS(:,end)',1.,'k-','LineWidth',0.2); hold on;
                    m_contour(lon_rho,lat_rho,mask_rho,[0.5 0.5],'k'); hold on;
                    m_grid('box','fancy','tickdir','in');
    figure(3); clf; plot(u_along(:,end),'r-'); hold on;
                    plot(v_normal(:,end),'b-'); hold on; legend('u along','v normal');
                    line([1:length(X)],repmat(0,1,length(X)),'Color','k');
    figure(4); clf; plot(u_along(:,end),'r-'); hold on;
                    plot(v_normal(:,end),'b-'); hold on; legend('u along','v normal');
                    line([1:length(X)],zeros(1,length(X)),'Color','k');
    figure(5); clf; plot(u_along(:,end),'r-'); hold on;
                    plot(v_normal(:,end),'b-'); hold on; legend('u along','v normal');
                    line([1:length(X)],repmat(0,1,length(X)),'Color','k');                
    clear lon_rho lat_rho mask_rho u v
  end  % if 0

  % Clean up
  clear rfil uin uin_EW vin vin_NS u_along

  % Set 0 currents on land and change variable name
  for i=1:length(X)
    for k=1:N
      if mask(i) < 0.5 | isnan(v_normal(i,k))
        u(i,k) = 0.;
      else
        u(i,k) = v_normal(i,k);
      end
    end
  end
  clear v_normal

  % Store u for use in calculating average along transect
  u_transect(:,:) = u;
  size(u)
  temp_transect(:,:) = temp;
  salt_transect(:,:) = salt;
  %no3_transect(:,:,fn) = no3_transect(:,:,fn) + no3;
  %phyto_transect(:,:,fn) = phyto_transect(:,:,fn) + phyto;
  disp(['Finished reading for ',int2str(dato(3)),'/',int2str(dato(2)),'-',int2str(dato(1))]);
  % Clean up
  clear dato tid u %Vin Vou
  % Find x-axis (sum of dx's)
  clear dxtot;
  dxtot(1) = 0.;
  for i=2:length(dx)
     dxtot(i) = dxtot(i-1) + dx(i);
  end
  dxtot = repmat(dxtot',1,N);

  % Define colors on colorbar
  b = [32 96 255;32 159 255;32 191 255;0 207 255;42 255 255;85 255 255;127 255 255;170 255 255; ...
       255 255 84;255 240 0;255 191 0;255 168 0;255 138 0;255 112 0;255 77 0;255 0 0]/255;

  % Define colorbar limits and color spacing
  mincol = -6.; maxcol = 6.; % Vlim = [-100.:(maxcol-mincol)/(length(b)-1):100.];

  %--------------------
  % Convert to Absolute Salinity and Conservative Temperature, 
  % calculate hydrostatic pressure and density.

  % Get pressure from water depth (negative in ocean)
  P = gsw_p_from_z( -z0, lat(1) );

  % Absolute Saliity
  SA = gsw_SA_from_SP( salt_transect, P, lon(1), lat(1) );

  % Conservative Temperature
  CT = gsw_CT_from_pt( SA, temp_transect );

  % Density
  RHO = gsw_rho( SA, CT, P );

  %--------------------

  % Fill contourf for transect :
  %   move bottom coordinates to sea bed and top coordinates to 0-level  
  z0h = [z0_w(:,1) z0(:,2:N-1) z0_w(:,N+1)];
  v_rho = 1026.5:0.5:1029;
  cd(rdir);
  % Plot umean1
  figure(1)
  clf
  %colormap(b)
  colormap(cm_balance); Vlim = [-100.:(maxcol-mincol)/(length(cm_balance)-1):100.];
  size(dxtot)
  [cs,hs] = contourf(dxtot/1000,-z0h,100.*u_transect,Vlim, ...
    'LineWidth',0.05, 'LineStyle','none');
  caxis([mincol maxcol]); 
  hold on;
  [cs2,hs2] = contour(dxtot/1000,-z0,RHO,v_rho, ...
    'LineColor','k','LineWidth',0.5,  ...
    'LineStyle','--');
  clabel(cs2,hs2,'FontSize',12,'LabelSpacing',300);
  cb = colorbar;
  cbarrow;
  %cb.Ticks = (mincol:(maxcol-mincol)/(length(b)-1):maxcol);
  cb.Label.String = 'cm/s';
  cb.Label.FontSize = 14;
  %cb_ytick(cb,mincol:(maxcol-mincol)/(length(b)-1):maxcol,1);
  %cb_title(cb,'cm/s',12);
  plot(dxtot(:,1)/1000,-h,'k-','LineWidth',2); hold on;
  text(0.03*dxtot(end,1)/1000 , 0.95*min(-h), 'South' , 'FontSize',14,'FontWeight','bold');
  text(0.85*dxtot(end,1)/1000 , 0.95*min(-h), 'North' , 'FontSize',14,'FontWeight','bold');
  xlabel('Distance (km)','FontSize',14,'FontWeight','bold');
  ylabel('Depth (m)','FontSize',14,'FontWeight','bold');
  set(gca,'FontSize',14,'FontWeight','bold');
  %writeVideo(vid, getframe(gcf));
  %%print('-dpng','-r600',['mean_curr_speed_JULAUGSEP_2005_WITH-FRESHWATER.png'])
  %print('-dpng','-r300',['mean_curr_speed_Kronebreen_K160v2_subglacial_JJAS_2007.png'])

  FigureName = strcat('curr_speed','_',MapFile,'_',num2str(aindex+1),'_',num2str(fn));
  hgsave(FigureName);
  %--------------------

  Vt = 0.0:0.1:6.0;

  % Plot Temperature
  figure(2)
  clf
  %colormap(b)
  [cs,hs] = contourf(dxtot/1000,-z0h,CT,Vt, ...
    'LineWidth',0.05, 'LineStyle','none');
  colormap(cm_thermal);
  %caxis([mincol maxcol]);
  caxis([3.0 5.0]);
  hold on;
  [cs2,hs2] = contour(dxtot/1000,-z0,RHO,v_rho, ...
    'LineColor','k','LineWidth',0.5,  ...
    'LineStyle','--');
  clabel(cs2,hs2,'FontSize',12,'LabelSpacing',300);
  cb = colorbar;
  cbarrow;
  %cb.Ticks = (mincol:(maxcol-mincol)/(length(b)-1):maxcol);
  cb.Label.String = '\Theta (^{\circ}C)';
  cb.Label.FontSize = 14;
  %cb_ytick(cb,mincol:(maxcol-mincol)/(length(b)-1):maxcol,1);
  %cb_title(cb,'cm/s',12);
  plot(dxtot(:,1)/1000,-h,'k-','LineWidth',2); hold on;
  text(0.03*dxtot(end,1)/1000 , 0.95*min(-h), 'South' , 'FontSize',14,'FontWeight','bold');
  text(0.85*dxtot(end,1)/1000 , 0.95*min(-h), 'North' , 'FontSize',14,'FontWeight','bold');
  xlabel('Distance (km)','FontSize',14,'FontWeight','bold');
  ylabel('Depth (m)','FontSize',14,'FontWeight','bold');
  set(gca,'FontSize',14,'FontWeight','bold');
  %print('-dpng','-r300','mean_temp_Kronebreen_K160v2_subglacial_JJAS_2007.png')
  FigureName = strcat('temp','_',MapFile,'_',num2str(aindex+1),'_',num2str(fn));
  hgsave(FigureName);

  %--------------------

  Vs = 19:0.1:35;
 
  % Plot Salinity
  figure(3)
  clf
  %colormap(b)
  [cs,hs] = contourf(dxtot/1000,-z0h,SA,Vs, ...
    'LineWidth',0.05, 'LineStyle','none');
  colormap(cm_haline);
  %caxis([mincol maxcol]);
  caxis([33 35]);
  hold on;
  [cs2,hs2] = contour(dxtot/1000,-z0,RHO,v_rho, ...
    'LineColor','k','LineWidth',0.5,  ...
    'LineStyle','--');
  clabel(cs2,hs2,'FontSize',12,'LabelSpacing',300);
  cb = colorbar;
  cbarrow;
  %cb.Ticks = (mincol:(maxcol-mincol)/(length(b)-1):maxcol);
  cb.Label.String = 'S_{A} (g / kg)';
  cb.Label.FontSize = 14;
  %cb_ytick(cb,mincol:(maxcol-mincol)/(length(b)-1):maxcol,1);
  %cb_title(cb,'cm/s',12);
  plot(dxtot(:,1)/1000,-h,'k-','LineWidth',2); hold on;
  text(0.03*dxtot(end,1)/1000 , 0.95*min(-h), 'South' , 'FontSize',14,'FontWeight','bold');
  text(0.85*dxtot(end,1)/1000 , 0.95*min(-h), 'North' , 'FontSize',14,'FontWeight','bold');
  xlabel('Distance (km)','FontSize',14,'FontWeight','bold');
  ylabel('Depth (m)','FontSize',14,'FontWeight','bold');
  set(gca,'FontSize',14,'FontWeight','bold');
  %print('-dpng','-r300','mean_salt_Kronebreen_K160v2_subglacial_JJAS_2007.png')
  FigureName = strcat('salt','_',MapFile,'_',num2str(aindex+1),'_',num2str(fn));
  hgsave(FigureName);
end  % fn
clear first fn Vstretching Vtransform theta_b theta_s hc 


