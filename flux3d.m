function [p2,p3,den1,info] = flux3d(dir0,nn,H,zr,zw)
% calculates terms in the 3D APE energy budget at a single time step
% then returns 3D fields [W m-3]


%% set values
% constants
g = 9.81;
offset = 1000; % to add to the density anomaly

% terms in the diagnostic balances
% rate or accel = sum of other terms
ts_list = {'hdiff','vdiff'};

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
f_avg = [dir0.avg,'ocean_avg_',ns,'.nc'];
f_dia = [dir0.dia,'ocean_dia_',ns,'.nc'];

% space and time info
%
% at average time (_avg)
[G,S,T] = Z_get_basic_info(f_avg);
ts_avg = T.ocean_time;
td_avg = T.time_datenum;
DZ = diff(zw);

%% Gather terms

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
[~,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,zr);
rho = den1 + offset; % use potential density
[D] = Z_flat(rho,rho,zr,zw,H);
zz = D.Z - D.Zf;


% form vertically-integrated PE terms
for ii = 1:length(ts_list)
    varname = ts_list{ii};
    p3.(varname) = g*zz.*rho.*DZ.* ...
        (-alpha.*t3.(varname) + beta.*s3.(varname));
    p2.(varname) = squeeze(sum(p3.(varname)));
end
% terms in p2 are [W m-2]


%% tidying up

info.ts_avg = ts_avg; info.td_avg = td_avg;

fnm = fieldnames(p2);
for ii = 1:length(fnm)
    varname = fnm{ii};
    p2.(varname)(H.island) = NaN;
end

