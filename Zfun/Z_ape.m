function [G, apea, apea_up, apea_down, zz, rr] = Z_ape(fn)
% Parker MacCready
%
% code to to calculate APE per unit area

[G,S,T] = Z_get_basic_info(fn);

island = Z_island(G);

[H] = Z_hyp(G,S,island);

%% calculate APE

% constants
g = 9.81;
rho0 = 1023.7; % set in the ROMS parameters
offset = 1000; % to add to the density anomaly

% space info
eta = nc_varget(fn,'zeta');
[z,z_w] = Z_s2z(G.h,eta,S);
DZ = diff(z_w);

%% intPE [W m-2]
salt = nc_varget(fn,'salt');
temp = nc_varget(fn,'temp');

[den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z);
rho = den1 + offset;

[D] = Z_flat_fast(rho,z,z_w,H);

zz = D.Z - D.Zf;
rr = D.R - D.Rf;

% divide the APE into up and down parts
zz_up = zz;
zz_down = zz;
zz_up(zz < 0) = 0;
zz_down(zz >= 0) = 0;

% create exact APE (neglecting compressibility)
F = g*(D.intRf - D.intRfzf);
apev = g*zz.*rho - F;
apev(apev < 0) = 0;
apea = squeeze(sum(apev.*DZ)); % [J m-2]

F_up = F;
F_up(zz < 0) = 0;
apev_up = g*zz_up.*rho - F_up;
apev_up(apev_up < 0) = 0;
apea_up = squeeze(sum(apev_up.*DZ)); % [J m-2]

F_down = F;
F_down(zz >= 0) = 0;
apev_down = g*zz_down.*rho - F_down;
apev_down(apev_down < 0) = 0;
apea_down = squeeze(sum(apev_down.*DZ)); % [J m-2]