close all; clear all; rehash;

% Load transect data
load('VolumeFlux_sub.mat');
load('VolumeFlux_surf.mat');
load('VolumeFlux_plm.mat');

DisEnt_KoB = csvread('Kongsbreen_output_86m\kongsbreen_results_disch_ent_n1.csv',1,0);
Discharge = squeeze(DisEnt_KoB(:,5));
Entrainment = squeeze(DisEnt_KoB(:,9));

figure
scatter(Discharge,Entrainment,8.0,'o');
hold on;
scatter(Discharge,(vf_KoB_sub - vf_KoB_surf)/10,5.0,'o');
%scatter(Discharge,Vflux_surf,5.0,'o');
%scatter(Discharge,Vflux_sub,5.0,'o');
legend('pyplume','ROMS / 10','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsbreen');

figure
scatter(Discharge,Entrainment,8.0,'o');
hold on;
scatter(Discharge(1:121),(vf_KoB_plm - vf_KoB_surf(1:121))/10,5.0,'o');
legend('pyplume','ROMS iceplume / 10','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsbreen');

figure
scatter(Discharge,vf_KoB_surf,8.0,'o');
hold on;
scatter(Discharge,vf_KoB_sub,8.0,'o');
set(gca,'yscale','log','YGrid','on');
legend('PSD','PTG','Location','SouthEast');
xlabel('Discharge volume flux');
ylabel('Inflow volume flux');
title('Kongsbreen');


DisEnt_KrB = csvread('Kronebreen_output\Kronebreen_results_disch_ent_nb.csv',1,0);
Discharge = squeeze(DisEnt_KrB(:,5));
Entrainment = squeeze(DisEnt_KrB(:,9));

figure
scatter(Discharge,Entrainment,8.0,'o');
hold on;
scatter(Discharge,(vf_KrB_sub - vf_KrB_surf)/13.5,5.0,'o');
%scatter(Discharge,Vflux_surf,5.0,'o');
%scatter(Discharge,Vflux_sub,5.0,'o');
legend('pyplume','ROMS / 13.5','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kronebreen');

figure
scatter(Discharge,Entrainment,8.0,'o');
hold on;
scatter(Discharge(1:121),(vf_KrB_plm - vf_KrB_surf(1:121))/13.5,5.0,'o');
legend('pyplume','ROMS iceplume / 13.5','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kronebreen');

DisEnt_KoV = csvread('Kongsvegen_output\Kongsvegen_results_disch_ent_nb.csv',1,0);
Discharge = squeeze(DisEnt_KoV(:,5));
Entrainment = squeeze(DisEnt_KoV(:,9));

figure
scatter(Discharge,Entrainment,8.0,'o');
hold on;
scatter(Discharge,(vf_KoV_sub - vf_KoV_surf)/17,5.0,'o');
%scatter(Discharge,Vflux_surf,5.0,'o');
%scatter(Discharge,Vflux_sub,5.0,'o');
legend('pyplume','ROMS / 17','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsvegen');

figure
scatter(Discharge,Entrainment,8.0,'o');
hold on;
scatter(Discharge(1:121),(vf_KoV_plm - vf_KoV_surf(1:121))/17,5.0,'o');
legend('pyplume','ROMS iceplume / 17','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsvegen');
