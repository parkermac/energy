function [apea, apea_up, apea_down] = Z_ape(fn,G,S,H)
% Parker MacCready
%
% Code to to calculate APE per unit area.  This is just the part associated
% with the density field deformations, and so we assume eta = 0.  The free
% surface contribution is added in later through the SW term.


%% calculate APE

% constants
g = 9.81;
offset = 1000; % to add to the density anomaly

% space info
eta = nc_varget(fn,'zeta');
[zr,zw] = Z_s2z(G.h,0*eta,S);
DZ = diff(zw);

% Debugging: print eta_flat
if 0
    island = Z_island(G);
    da = G.DX .* G.DY;
    eta_flat = Z_eta_flat(eta, da, island);
    format long
    disp(['eta_flat = ',sprintf('%.12f',eta_flat)])
end

%% intPE [W m-2]
salt = nc_varget(fn,'salt');
temp = nc_varget(fn,'temp');

[~,den1,~,~] = Z_roms_eos_vec(salt,temp,zr);
rho = den1 + offset; % use potential density

[D] = Z_flat(rho,rho,zr,zw,H);

zz = D.Z - D.Zf;

% divide the APE into up and down parts
zz_up = zz;
zz_down = zz;
zz_up(zz < 0) = 0;
zz_down(zz >= 0) = 0;

% create exact APE (neglecting compressibility)
F = g*(D.intRf - D.intRfzf);
apev = g*zz.*rho - F;
if 0
    disp('Warning: making apev non-negative')
    disp(['  Number of negative points = ',num2str(nansum(apev(:) < 0))])
    disp(['  Mean of negative points = ',num2str(nanmean(apev((apev(:) < 0))))])
    disp(['  Largest of negative points = ',num2str(nanmin(apev((apev(:) < 0))))])
    disp(['  Mean of all points = ',num2str(nanmean(apev(:)))])
end
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