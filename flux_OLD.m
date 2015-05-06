function [p2,k2,info] = flux_OLD(dir0,nn,H)
% 6/12/2014  Parker MacCready
%
% calculates terms in the 3D energy budget at a single time step
% then returns vertical integrals [W m-2]
%
% * carefully zeros-out momentum terms to retain KE momentum terms over
%   the full rho-grid
%
% * uses the exact ROMS compressible equation of state, although we use
%   the density referenced to surface pressure for all APE calculations
%
% * calculates an approximation of the APE terms from
%   Scotti & White (2014).  This requires us to specify "island"
%   inside the code.
%
% * uses hdiff instead of xdiff + ydiff for the tracer advection, which
%   leads to nearly exact APE balance (see dia_DEBUG.m)

%% set values

% constants
g = 9.81;
g2 = g/2;
rho0 = 1023.7; % set in the ROMS parameters
offset = 1000; % to add to the density anomaly

do_scotti = 1; % a flag passed to Z_flat.m

% terms in the diagnostic balances
% rate or accel = sum of other terms
ts_list = {'rate','xadv','yadv','vadv','hdiff','vdiff'};
uv_list = {'accel','xadv','yadv','vadv','cor','prsgrd','vvisc'};

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
f_his1 = [dir0.his,'ocean_his_',ns,'.nc'];
f_avg = [dir0.avg,'ocean_avg_',ns,'.nc'];
f_dia = [dir0.dia,'ocean_dia_',ns,'.nc'];
ns2 = num2str(nn+1); ns2 = ['0000',ns2]; ns2 = ns2(end-3:end);
f_his2 = [dir0.his,'ocean_his_',ns2,'.nc'];
% Note: his1 and his2 bracket the averaging period for avg & dia

% space and time info
%
% at average time (_avg)
[G,S,T] = Z_get_basic_info(f_avg);
ts_avg = T.ocean_time;
td_avg = T.time_datenum;
eta_avg = nc_varget(f_avg,'zeta');
[z_avg,z_w_avg] = Z_s2z(G.h,eta_avg,S);
DZ_avg = diff(z_w_avg);
% at start and end times (1, 2)
T = Z_get_time_info(f_his1);
t1 = T.ocean_time; % time in seconds
T = Z_get_time_info(f_his2);
t2 = T.ocean_time; % time in seconds
DT = t2 - t1;

%% Potential Energy calculation

% get 3D tracer diagnostics, and save into structures
% temperature [C s-1] & salinity [psu s-1] (rho grid)
for ii = 1:length(ts_list)
    varname = ts_list{ii};
    tname = ['temp_',varname];
    sname = ['salt_',varname];
    t3.(varname) = nc_varget(f_dia,tname);
    s3.(varname) = nc_varget(f_dia,sname);
end

salt = nc_varget(f_avg,'salt');
temp = nc_varget(f_avg,'temp');

% see written notes for explanation
[den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z_avg);
rho_avg = den1 + offset;
rho00 = rho_avg.*(1 + alpha.*temp - beta.*salt);
[D_avg] = Z_flat(eta_avg,rho_avg,z_avg,z_w_avg,H,do_scotti);
zz_avg = D_avg.Z - D_avg.Zf;
rr_avg = D_avg.R - D_avg.Rf;
clear salt temp

% form vertically-integrated PE terms
for ii = 1:length(ts_list)
    varname = ts_list{ii};
    p3.(varname) = -alpha.*t3.(varname) + beta.*s3.(varname);
    p2.(varname) = squeeze(sum(g2*zz_avg.*rho_avg.*DZ_avg.*p3.(varname)));
end
% terms in p2 are [W m-2]

clear t3 s3 p3

% also save the APE itself, since we have all the pieces
p2.ape = squeeze(sum(g2*zz_avg.*DZ_avg.*rr_avg)); % [J m-2]

% we have to add a term to the rate to get it to represent the
% quantity we want, and this extra term will be added to the
% advective terms on the RHS to ensure balance

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

% then make additional PE terms to make storage rate meaningful
%
dpedt1 = squeeze(sum(g2*zz_avg.*rho00.*(DZ2 - DZ1)/DT));
dpedt2 = squeeze(sum(-g2*zz_avg.*(D2.Rf.*DZ2 - D1.Rf.*DZ1)/DT));
dpedt3 = squeeze(sum(g2*rr_avg.*DZ_avg.*(zz2 - zz1)/DT));

% and add them to the results
p2.rate = p2.rate + dpedt1 + dpedt2 + dpedt3;
p2.vadv = p2.vadv + dpedt1;
p2.dpe0dt = dpedt2 + dpedt3;

% also create dAPE/dt from scratch to check our calculation
%
pe1 = squeeze(sum(g2*zz1.*DZ1.*rr1));
pe2 = squeeze(sum(g2*zz2.*DZ2.*rr2));
p2.rate_check = (pe2 - pe1)/DT;

% cleaning up
clear pe1 pe2
clear dpedt1 dpedt2 dpedt3
clear rho_avg rho1 rho2
clear zz_avg zz1 zz2 rr_avg rr1 rr2 D_avg D1 D2
clear alpha beta den den1

%% Kinetic Energy calculation

% get 3D momentum diagnostics, and save into structures [m s-2]
% (u grid and v grid)
for ii = 1:length(uv_list)
    varname = uv_list{ii};
    uname = ['u_',varname];
    vname = ['v_',varname];
    u3.(varname) = nc_varget(f_dia,uname);
    v3.(varname) = nc_varget(f_dia,vname);
end

% form vertically-integrated KE terms [W m-2]
% (still u and v grids)
Huon = nc_varget(f_avg,'Huon');
Hvom = nc_varget(f_avg,'Hvom');
DY_u = interp2(G.lon_rho,G.lat_rho,G.DY,G.lon_u,G.lat_u);
DX_v = interp2(G.lon_rho,G.lat_rho,G.DX,G.lon_v,G.lat_v);
for ii = 1:length(uv_list)
    varname = uv_list{ii};
    ku3.(varname) = rho0*Huon.*u3.(varname);
    kv3.(varname) = rho0*Hvom.*v3.(varname);
    ku2.(varname) = squeeze(sum(ku3.(varname)))./DY_u;
    kv2.(varname) = squeeze(sum(kv3.(varname)))./DX_v;
end

clear u3 v3 ku3 kv3

% create and apply fix terms
u1 = nc_varget(f_his1,'u');
u2 = nc_varget(f_his2,'u');
v1 = nc_varget(f_his1,'v');
v2 = nc_varget(f_his2,'v');
ut = (u2 - u1)/DT;
vt = (v2 - v1)/DT;
accel_fix_u = -0.5*ku2.accel + 0.5*rho0*squeeze(sum(Huon.*ut))./DY_u;
accel_fix_v = -0.5*kv2.accel + 0.5*rho0*squeeze(sum(Hvom.*vt))./DX_v;
% add fix terms to LHS and RHS
ku2.accel = ku2.accel + accel_fix_u;
kv2.accel = kv2.accel + accel_fix_v;
ku2.xadv = ku2.xadv + accel_fix_u;
kv2.yadv = kv2.yadv + accel_fix_v;

clear u1 u2 v1 v2 ut vt
clear Huon Hvom DY_u DX_v
clear accel_fix_u accel_fix_v

% interpolate to rho grid
for ii = 1:length(uv_list)
    varname = uv_list{ii};
    % note that we set values to zero on land
    % before interpolating to rho grid
    ku2.(varname)(~G.mask_u) = 0;
    kv2.(varname)(~G.mask_v) = 0;
    k2.(varname) = ...
        interp2(G.lon_u,G.lat_u,ku2.(varname),G.lon_rho,G.lat_rho) + ...
        interp2(G.lon_v,G.lat_v,kv2.(varname),G.lon_rho,G.lat_rho);
end

% form d(KE)/dt from scratch [W m-2]
ke1 = Z_make_ke(G,S,rho0,f_his1);
ke2 = Z_make_ke(G,S,rho0,f_his2);
k2.accel_check = (ke2 - ke1)/DT;
clear ke2 ke1

% form the KE reservoir [J m-2]
k2.ke = Z_make_ke(G,S,rho0,f_avg);

%% tidying up

% save info
ts_list{length(ts_list)+1} = 'dpe0dt';
ts_list{length(ts_list)+1} = 'rate_check';
ts_list{length(ts_list)+1} = 'ape';
uv_list{length(uv_list)+1} = 'accel_check';
uv_list{length(uv_list)+1} = 'ke';
info.ts_list = ts_list; info.uv_list = uv_list;
info.ts_avg = ts_avg; info.td_avg = td_avg;

for vv = 1:length(uv_list)
    varname = uv_list{vv};
    k2.(varname)(H.island) = NaN;
end
for vv = 1:length(ts_list)
    varname = ts_list{vv};
    p2.(varname)(H.island) = NaN;
end

