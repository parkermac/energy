function [ex2,info] = extras(dir0,nn,H,G,S)
% calculates extra terms in the 3D energy budget at a single time step
% then returns vertical integrals [W m-2]

%% set values

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
f_avg = [dir0.avg,'ocean_avg_',ns,'.nc'];
f_dia = [dir0.dia,'ocean_dia_',ns,'.nc'];

% space and time info
%
% at average time (_avg)
T = Z_get_time_info(f_avg);
ts_avg = T.ocean_time;
td_avg = T.time_datenum;

%% Wind work on geostrophic flow

% pressure gradients
% assumed to be of the form, e.g.: -dpdx/rho0 [m s-2]
px = nc_varget(f_dia,'u_prsgrd',[0, S.N-1, 0, 0],[1, 1, -1, -1]);
py = nc_varget(f_dia,'v_prsgrd',[0, S.N-1, 0, 0],[1, 1, -1, -1]);

% geostrophic velicity
vg = -px./sw_f(G.lat_u);
ug = py./sw_f(G.lat_v);

VG = interp2(G.lon_u,G.lat_u,vg,G.lon_v,G.lat_v);
UG = interp2(G.lon_v,G.lat_v,ug,G.lon_u,G.lat_u);

% form the wind work [W m-2]
sustr = squeeze(nc_varget(f_avg,'sustr',[0, 0, 0],[1, -1, -1]));
svstr = squeeze(nc_varget(f_avg,'svstr',[0, 0, 0],[1, -1, -1]));
ex2.wind_work =  ...
        interp2(G.lon_u,G.lat_u,UG.*sustr,G.lon_rho,G.lat_rho) + ...
        interp2(G.lon_v,G.lat_v,VG.*svstr,G.lon_rho,G.lat_rho);

%% tidying up

info.ts_avg = ts_avg; info.td_avg = td_avg;

fnm = fieldnames(ex2);
for ii = 1:length(fnm)
    varname = fnm{ii};
    ex2.(varname)(H.island) = NaN;
end
