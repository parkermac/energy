% Code to plot the density and the flattened version.
% First developed for WTD 2015.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');
fn = '/Users/PM5/Documents/roms/output/D2005_avg/ocean_avg_1836.nc';

%% get info
[G,S,T] = Z_get_basic_info(fn);
island = Z_island(G);
[H] = Z_hyp(G,S,island);

% constants
g = 9.81;
rho0 = 1023.7; % set in the ROMS parameters
offset = 1000; % to add to the density anomaly

% space info
eta = nc_varget(fn,'zeta');
[z,z_w] = Z_s2z(G.h,eta,S);
DZ = diff(z_w);

% density
salt = nc_varget(fn,'salt');
temp = nc_varget(fn,'temp');
[den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z);
rho = den1 + offset;

% flattened state
[D,Zf,Rf] = Z_flat_fast(rho,z,z_w,H);

zz = D.Z - D.Zf;
rr = D.R - D.Rf;

% divide the APE into up and down parts
zz_up = zz;
zz_down = zz;
zz_up(zz < 0) = 0;
zz_down(zz >= 0) = 0;

%% plotting
close all
Z_fig(18)
set(gcf,'position',[10 10 1500 750]);

ny = 200;
x = cumsum(G.DX(ny,:))/1000;
x_sect_full = ones(S.N + 2,1) * x(:)';
z_sect = squeeze(D.Z(:,ny,:));
z_sect_full = [-G.h(ny,:); z_sect; eta(ny,:)];
sigma_sect = squeeze(D.R(:,ny,:)) - 1000;
sigma_sect_full = [sigma_sect(1,:); sigma_sect; sigma_sect(end,:)];
sigmaf_sect = squeeze(D.Rf(:,ny,:)) - 1000;
sigmaf_sect_full = [sigmaf_sect(1,:); sigmaf_sect; sigmaf_sect(end,:)];

sig_lims = [23, 27];
x_lims = [40, 250];
z_lims = [-500, 5];

subplot(131)
pcolor(x_sect_full,z_sect_full,sigma_sect_full);
shading interp
caxis(sig_lims);
colormap jet
%colorbar('south')
axis([x_lims, z_lims]);
title('\sigma_{0} (kg m^{-3})');
xlabel('Zonal Distance (km)');
ylabel('Z (m)');

subplot(132)
pcolor(x_sect_full,z_sect_full,sigmaf_sect_full);
shading interp
caxis(sig_lims);
colormap jet
%colorbar('south')
set(gca,'YTickLabel',[])
axis([x_lims, z_lims]);
title('Flattened \sigma_{0} (kg m^{-3})');
xlabel('Zonal Distance (km)');

subplot(133)
lh = plot(D.R(:) - 1000,D.Z(:),'.c',Rf - 1000,Zf,'.b','markersize',1);
set(gca,'YTickLabel',[])
axis([sig_lims, z_lims]);
title('Density Profiles')
xlabel('\sigma_{0} (kg m^{-3})');
[xt,yt] = Z_lab('ll');
text(xt,yt + 30,'All Points in Volume','color','c')
text(xt,yt,'Flattened State','color','b')

%%
if 1
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',[Tdir.output,'energy_out/density_plot.png']);
end


