% APE_test.m  6/12/2014  Parker MacCready
%
% code to test different methods to calculate APE per unit area
%
% Added lines to make apev positive definite.  There were a few very small
% negative points in the test case (min was -1.8e-9 compared to a median of
% +6 and a mean of +77 in apev).  This could be due to roundoff error.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

basename = 'D2005';
nn = 1836; %
dir1 = '/Users/PM5/Documents/roms/output/';
dir0.avg = [dir1,basename,'_avg/'];
dir0.dia = [dir1,basename,'_dia/'];
dir0.his = [dir1,basename,'_his/'];
%% get hypsometry structure H

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
%fn = [dir0.his,'ocean_his_',ns,'.nc'];
fn = [dir0.avg,'ocean_avg_',ns,'.nc'];
[G,S,T] = Z_get_basic_info(fn);

island = Z_island(G);

[H] = Z_hyp(G,S,island);

%% calculate APE

% constants
g = 9.81;
rho0 = 1023.7; % set in the ROMS parameters
offset = 1000; % to add to the density anomaly

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
%rho00 = rho_avg.*(1 + alpha.*temp - beta.*salt);
[D_avg] = Z_flat(eta_avg,rho_avg,rho_avg,z_avg,z_w_avg,H);
zz_avg = D_avg.Z - D_avg.Zf;
rr_avg = D_avg.R - D_avg.Rf;

% divide the APE into up and down parts
zz_up = zz_avg;
zz_down = zz_avg;
zz_up(zz_avg < 0) = 0;
zz_down(zz_avg >= 0) = 0;

% create exact APE (neglecting compressibility)
F_avg = g*(D_avg.intRf - D_avg.intRfzf);
apev = g*zz_avg.*rho_avg - F_avg;
apev(apev < 0) = 0;
apea = squeeze(sum(apev.*DZ_avg)); % [J m-2]

F_up = F_avg;
F_up(zz_avg < 0) = 0;
apev_up = g*zz_up.*rho_avg - F_up;
apev_up(apev_up < 0) = 0;
apea_up = squeeze(sum(apev_up.*DZ_avg)); % [J m-2]

F_down = F_avg;
F_down(zz_avg >= 0) = 0;
apev_down = g*zz_down.*rho_avg - F_down;
apev_down(apev_down < 0) = 0;
apea_down = squeeze(sum(apev_down.*DZ_avg)); % [J m-2]

% test of the split
apea_alt = apea - apea_up - apea_down;
% RESULT: this is essentially zero, so the split is good


%% plotting
close all
figure
Z_fig(14)
set(gcf,'position',[10 10 2500 1000]);
cax = [2 6];


fld_list = {'apea','apea_up','apea_down'};

for ii = 1:length(fld_list)
    
    fld_name = fld_list{ii};
    
    eval(['fld = ',fld_name,';']);
    
    fld(fld<0) = nan; % there is one bad point in the SW corner...?
    
    subplot(1,3,ii)
    colormap jet
    Z_pcolorcen(G.lon_rho,G.lat_rho,log10(fld));
    caxis(cax);
    shading flat
    colorbar
    hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
    title(['log10 ',strrep(fld_name,'_',' '),' (J m^{-2})'])
    Z_dar;
    Z_addcoast('combined',Tdir.coast);
    [xt,yt] = Z_lab('ll');
    text(xt,yt,datestr(T.time_datenum))

end


