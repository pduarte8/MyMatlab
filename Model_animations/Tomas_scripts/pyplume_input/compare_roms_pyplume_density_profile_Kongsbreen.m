close all; clear all; rehash;

% Load plume data
Zpp = csvread('Kongsbreen_output_86m\kongsbreen_results_z.csv');
Tpp = csvread('Kongsbreen_output_86m\kongsbreen_results_t_p.csv');
Spp = csvread('Kongsbreen_output_86m\kongsbreen_results_s_p.csv');
Ta = csvread('Kongsbreen_output_86m\kongsbreen_results_t_a.csv');
Sa = csvread('Kongsbreen_output_86m\kongsbreen_results_s_a.csv');

Hpp = squeeze(Zpp(1,:) - Zpp(1,end));

Ppp = gsw_p_from_z( Hpp, 79.0 );

RHOpp = gsw_rho( Spp, Tpp, Ppp );
RHOa = gsw_rho( Sa, Ta, Ppp );

load('Kongsbreen_profile.mat');

load('Kongsbreen_amb_profile.mat');

for i=1:122
  %plot(RHO_KoB(i,:),z_KoB,'b-',RHOpp(i,:),-Zpp(i,:),'r-',RHOa(i,:),-Zpp(i,:),'g-');
  plot(RHO_KoB(i,:),z_KoB,'b-',RHO_KoB_amb(i,:),z_KoB_amb,'g-',RHOpp(i,:),Hpp,'r-');
  %plot(RHO_KoB(i,:),z_KoB,'b-',RHO_KoB_amb(i,:),z_KoB_amb,'g-');
  axis([1000 1029 -120 0]);
  %pause(0.2);
end

figure
p1 = gsw_SA_CT_plot(SA_KoB_amb,CT_KoB_amb);
p2 = gsw_SA_CT_plot(SA_KoB,CT_KoB);
set(p2,'Marker','o','MarkerSize',2,'MarkerFaceColor','r','MarkerEdgeColor','none');
p3 = gsw_SA_CT_plot(Spp,Tpp);
set(p3,'Marker','o','MarkerSize',2,'MarkerFaceColor','g','MarkerEdgeColor','none');

DisEnt_KoB = csvread('Kongsbreen_output_86m\kongsbreen_results_disch_ent_n1.csv',1,0);
[Discharge,I] = sort(squeeze(DisEnt_KoB(:,5)));
RHO_Ambient = RHO_KoB_amb(I,:);
RHO_roms = RHO_KoB(I,:);
RHO_pyplume = RHOpp(I,:);
RHO_pyp_amb = RHOa(I,:);
