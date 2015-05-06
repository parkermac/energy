% APE_test.m  6/12/2014  Parker MacCready
%
% code to test different methods to calculate APE per unit area

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

basename = 'D2005';
nn = 1848; %
dir1 = '/Users/PM5/Documents/roms/output/';
dir0.avg = [dir1,basename,'_avg/'];
dir0.dia = [dir1,basename,'_dia/'];
dir0.his = [dir1,basename,'_his/'];
%% get hypsometry structure H

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
fn = [dir0.his,'ocean_his_',ns,'.nc'];
[G,S,T] = Z_get_basic_info(fn);

island = Z_island(G);

[H] = Z_hyp(G,S,island);

%% calculate APE in different ways

% constants
g = 9.81;
g2 = g/2;
rho0 = 1023.7; % set in the ROMS parameters
offset = 1000; % to add to the density anomaly

do_scotti = 1;

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
f_his1 = [dir0.his,'ocean_his_',ns,'.nc'];
f_avg = [dir0.avg,'ocean_avg_',ns,'.nc'];
f_dia = [dir0.dia,'ocean_dia_',ns,'.nc'];
ns2 = num2str(nn+1); ns2 = ['0000',ns2]; ns2 = ns2(end-3:end);
f_his2 = [dir0.his,'ocean_his_',ns2,'.nc'];
% Note: his1 and his2 bracket the averaging period for avg & dia

% space info
[G,S,T] = Z_get_basic_info(f_avg);
ts_avg = T.ocean_time;
td_avg = T.time_datenum;
eta_avg = nc_varget(f_avg,'zeta');
[z_avg,z_w_avg] = Z_s2z(G.h,eta_avg,S);
DZ_avg = diff(z_w_avg);

%% intPE [W m-2]
salt = nc_varget(f_avg,'salt');
temp = nc_varget(f_avg,'temp');

[den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z_avg);
rho_avg = den1 + offset;
rho00 = rho_avg.*(1 + alpha.*temp - beta.*salt);
[D_avg] = Z_flat(eta_avg,rho_avg,z_avg,z_w_avg,H,do_scotti);
zz_avg = D_avg.Z - D_avg.Zf;
rr_avg = D_avg.R - D_avg.Rf;

% create the approximate APE
ape_avg = squeeze(sum(g2*zz_avg.*DZ_avg.*rr_avg)); % [J m-2]

% create exact APE (neglecting compressibility)
F_avg = g*(D_avg.intRf - D_avg.intRfzf);
apev = g*zz_avg.*rho_avg - F_avg;
ape_alt = squeeze(sum(apev.*DZ_avg)); % [J m-2]
ape_alt(ape_alt<=0) = 1e-10; % there is one negative point, with value -2.8e-7

%% testing
% first get some start and end values
%
eta1 = nc_varget(f_his1,'zeta');
[z1,z_w] = Z_s2z(G.h,eta1,S);
DZ1 = diff(z_w);
salt = nc_varget(f_his1,'salt');
temp = nc_varget(f_his1,'temp');
[den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z1);
rho1 = den1 + offset;
clear salt temp
[D1] = Z_flat(eta1,rho1,z1,z_w,H,do_scotti);
zz1 = D1.Z - D1.Zf;
rr1 = D1.R - D1.Rf;
%
eta2 = nc_varget(f_his2,'zeta');
[z2,z_w] = Z_s2z(G.h,eta2,S);
DZ2 = diff(z_w);
salt = nc_varget(f_his2,'salt');
temp = nc_varget(f_his2,'temp');
[den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z2);
rho2 = den1 + offset;
clear salt temp
[D2] = Z_flat(eta2,rho2,z2,z_w,H,do_scotti);
zz2 = D2.Z - D2.Zf;
rr2 = D2.R - D2.Rf;

% at start and end times (1, 2)
T = Z_get_time_info(f_his1);
t1 = T.ocean_time; % time in seconds
T = Z_get_time_info(f_his2);
t2 = T.ocean_time; % time in seconds
DT = t2 - t1;

drho_dt = (rho2-rho1)/DT;
[D_dot] = Z_flat(eta_avg,drho_dt,z_avg,z_w_avg,H,do_scotti);


%% plotting
close all
figure
Z_fig(14)
set(gcf,'position',[10 10 800 600]);
cax = [1 6];

subplot(131)
Z_pcolorcen(G.lon_rho,G.lat_rho,log10(ape_avg));
caxis(cax);
shading flat
colorbar
title('log10 APE (J m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);
[xt,yt] = Z_lab('ll');
text(xt,yt,datestr(T.time_datenum))

subplot(132)
Z_pcolorcen(G.lon_rho,G.lat_rho,log10(ape_alt));
caxis(cax);
shading flat
colorbar
title('log10 APE alt (J m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);
[xt,yt] = Z_lab('ll');

subplot(133)
Z_pcolorcen(G.lon_rho,G.lat_rho,log10(abs(ape_avg - ape_alt)));
caxis(cax);
shading flat
colorbar
title('log10 abs Difference (J m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);
[xt,yt] = Z_lab('ll');

