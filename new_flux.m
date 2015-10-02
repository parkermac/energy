function p2 = new_flux(dir0,nn,H)
% calculates dAPEa/dt at a single time step [W m-2]
%
% The purpose is to better understand the APE calculation.

%% set values

% constants
g = 9.81;
rho0 = 1023.7; % set in the ROMS parameters
offset = 1000; % to add to the density anomaly

% terms in the diagnostic balances
% rate or accel = sum of other terms
ts_list = {'rate'};

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
[z_avg,z_w_avg] = Z_s2z(G.h,0*eta_avg,S);
DZ_avg = diff(z_w_avg);
% at start and end times (1, 2)
T = Z_get_time_info(f_his1);
t1 = T.ocean_time; % time in seconds
T = Z_get_time_info(f_his2);
t2 = T.ocean_time; % time in seconds
DT = t2 - t1;

%% Potential Energy calculation
% Note: the APE calculated here does not include the contribution of the
% free surface.  This is done later in Z_make_derived.m by adding sw.ape.

% get 3D tracer diagnostics, and save into structures
% temperature [C s-1] & salinity [psu s-1] (rho grid)
for ii = 1:length(ts_list)
    varname = ts_list{ii};
    tname = ['temp_',varname];
    sname = ['salt_',varname];
    t3.(varname) = nc_varget(f_dia,tname);
    s3.(varname) = nc_varget(f_dia,sname);
end

salt_avg = nc_varget(f_avg,'salt');
temp_avg = nc_varget(f_avg,'temp');
% see written notes for explanation
[~,den1,alpha_avg,beta_avg] = Z_roms_eos_vec(salt_avg,temp_avg,z_avg);
rho_avg = den1 + offset;
% NOTE: the alpha and beta are
% in situ values, whereas the density we use is referenced to zero
% pressure.  The difference of these terms for s=31, theta=8 between 0 and
% 1000 m is 1.2% for beta, but 14% for alpha. 5/4/2015
% It is not clear which is correct.
[D_avg] = Z_flat(rho_avg,rho_avg,z_avg,z_w_avg,H);
zz_avg = D_avg.Z - D_avg.Zf;

% form vertically-integrated PE terms
for ii = 1:length(ts_list)
    varname = ts_list{ii};
    p3.(varname) = -alpha_avg.*t3.(varname) + beta_avg.*s3.(varname);
    p2.(varname) = squeeze(sum(g*zz_avg.*rho_avg.*DZ_avg.*p3.(varname)));
end
% terms in p2 are [W m-2]

clear t3 s3 p3

% We have to add some terms to get the LHS to be the rate that
% we want.  These will be added to the RHS to ensure balance.

% first get some start and end values
%
eta1 = nc_varget(f_his1,'zeta');
[z1,z_w] = Z_s2z(G.h,eta1,S);
DZ1 = diff(z_w);
salt = nc_varget(f_his1,'salt');
temp = nc_varget(f_his1,'temp');
[~,den1,~,~] = Z_roms_eos_vec(salt,temp,z1);
rho1 = den1 + offset;
clear salt temp
[D1] = Z_flat(rho1,rho1,z1,z_w,H);
%
eta2 = nc_varget(f_his2,'zeta');
[z2,z_w] = Z_s2z(G.h,eta2,S);
DZ2 = diff(z_w);
salt = nc_varget(f_his2,'salt');
temp = nc_varget(f_his2,'temp');
[~,den1,~,~] = Z_roms_eos_vec(salt,temp,z2);
rho2 = den1 + offset;
clear salt temp
[D2] = Z_flat(rho2,rho2,z2,z_w,H);

% then make additional PE terms to make storage rate meaningful

% rate stuff associated with delta_t
albet = alpha_avg.*temp_avg - beta_avg.*salt_avg;
p2.dapedt_uno = squeeze(sum(g*zz_avg.*rho_avg.*(albet).*(DZ2 - DZ1)/DT));

% term associated with rate of change of rho_flat.
fld = (D2.Rf - D1.Rf)/DT;
[D_dot] = Z_flat(rho_avg,fld,z_avg,z_w_avg,H);
F_dot = g*(D_dot.intFf - D_dot.intFfzf);
p2.dapedt_dos = squeeze(sum(-F_dot.*DZ_avg));

% and add them to the results
p2.rate = p2.rate + p2.dapedt_uno;
p2.rate = p2.rate + p2.dapedt_dos;

% also create dAPEa/dt from scratch to check our calculation
[apea1, ~, ~] = Z_ape(f_his1,G,S,H);
[apea2, ~, ~] = Z_ape(f_his2,G,S,H);
p2.rate_check = (apea2 - apea1)/DT;

% cleaning up
clear pe1 pe2
clear dpedt1 dpedt2 dpedt3
clear rho_avg rho1 rho2
clear zz_avg zz1 zz2 rr_avg rr1 rr2 D_avg D1 D2
clear alpha beta den1
clear apea1 apea2
clear F1 F2 F_dot dapedt_uno dapedt_dos albet

%% tidying up

fnm = fieldnames(p2);
for ii = 1:length(fnm)
    varname = fnm{ii};
    p2.(varname)(H.island) = NaN;
end

