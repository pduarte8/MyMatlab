close all; clear all; rehash;

% Load transect data
load('volflux_Kongsbreen/Transect_sub.mat');
load('volflux_Kongsbreen/Transect_surf.mat');

% Load plume data
Zpp = csvread('Kongsbreen_output_86m\kongsbreen_results_z.csv');
Tpp = csvread('Kongsbreen_output_86m\kongsbreen_results_t_p.csv');
Spp = csvread('Kongsbreen_output_86m\kongsbreen_results_s_p.csv');
Ta = csvread('Kongsbreen_output_86m\kongsbreen_results_t_a.csv');
Sa = csvread('Kongsbreen_output_86m\kongsbreen_results_s_a.csv');

DisEnt_KoB = csvread('Kongsbreen_output_86m\kongsbreen_results_disch_ent_n1.csv',1,0);
[Discharge,I] = sort(squeeze(DisEnt_KoB(:,5)));
Entrainment = squeeze(DisEnt_KoB(I,9));
Vflux_sub = Vin_sub(I);
Vflux_surf = Vin_surf(I);

figure
scatter(Discharge,Entrainment,8.0,'o');
hold on;
scatter(Discharge,(Vflux_sub - Vflux_surf)/12,5.0,'o');
%scatter(Discharge,Vflux_surf,5.0,'o');
%scatter(Discharge,Vflux_sub,5.0,'o');
legend('pyplume','ROMS / 12','Location','NorthWest');
xlabel('Discharge volume flux');
ylabel('Entrainment volume flux');
title('Kongsbreen');
