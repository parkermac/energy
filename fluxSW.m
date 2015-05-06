function [sw,sw_list] = fluxSW(dir0,nn,island)
% 6/12/2014  Parker MacCready
% calculates terms in the SW energy budget at a single time step [W m-2]
%
% set all masked areas (u, v grids) to 0 instead of NaN

%% set values of constants

g = 9.81;
rho0 = 1023.7; % set in the ROMS parameters

uv_list = {'accel','xadv','yadv','cor','prsgrd','sstr','bstr'};

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
f_his1 = [dir0.his,'ocean_his_',ns,'.nc'];
f_avg = [dir0.avg,'ocean_avg_',ns,'.nc'];
f_dia = [dir0.dia,'ocean_dia_',ns,'.nc'];
ns2 = num2str(nn+1); ns2 = ['0000',ns2]; ns2 = ns2(end-3:end);
f_his2 = [dir0.his,'ocean_his_',ns2,'.nc'];
% Note: his1 and his2 bracket the averaging period for avg & dia
[G,S,T] = Z_get_basic_info(f_avg);
eta_avg = nc_varget(f_avg,'zeta');

%% Energy calculation

% MASS terms
%
% form deta/dt (or dh/dt) from history files
e1 = nc_varget(f_his1,'zeta');
e2 = nc_varget(f_his2,'zeta');
%
t1 = nc_varget(f_his1,'ocean_time');
t2 = nc_varget(f_his2,'ocean_time');
DT = t2 - t1;
dedt = (e2 - e1)/DT; % = deta/dt [m s-1]
ee1 = e1.*e1; ee2 = e2.*e2;
deedt = (ee2 - ee1)/DT; % d(eta^2)/dt [m2 s-1]
%
% flattened surface heights
da = G.DX .* G.DY;
eaf = sum(eta_avg(~island).*da(~island))/sum(da(~island));
e1f = sum(e1(~island).*da(~island))/sum(da(~island));
e2f = sum(e2(~island).*da(~island))/sum(da(~island));
defdt = (e2f - e1f)/DT; % = def/dt [m s-1]
%
% differences from flattened state (p for prime)
epa = eta_avg - eaf;
ep1 = e1 - e1f;
ep2 = e2 - e2f;
depdt = (ep2 - ep1)/DT;
eepa = epa.*epa;
eep1 = ep1.*ep1;
eep2 = ep2.*ep2;
deepdt = (eep2 - eep1)/DT;

% also du/dt & dvdt
% u_avg = nc_varget(f_avg,'ubar');
% v_avg = nc_varget(f_avg,'vbar');
u1 = nc_varget(f_his1,'ubar');
u2 = nc_varget(f_his2,'ubar');
v1 = nc_varget(f_his1,'vbar');
v2 = nc_varget(f_his2,'vbar');
dudt = (u2 - u1)/DT;
dvdt = (v2 - v1)/DT;
dudt(~G.mask_u) = 0;
dvdt(~G.mask_v) = 0;

% form hu and hv from averages
Huon = nc_varget(f_avg,'Huon');
Hvom = nc_varget(f_avg,'Hvom');
hudy = squeeze(sum(Huon)); % [m3 s-1] (u grid)
hvdx = squeeze(sum(Hvom)); % [m3 s-1] (v grid)
hudy(~G.mask_u) = 0;
hvdx(~G.mask_v) = 0;
DY_u = interp2(G.lon_rho,G.lat_rho,G.DY,G.lon_u,G.lat_u);
DX_v = interp2(G.lon_rho,G.lat_rho,G.DX,G.lon_v,G.lat_v);
hu = hudy ./ DY_u; % [m2 s-1] (u grid)
hv = hvdx ./ DX_v; % [m2 s-1] (v grid)

% rate of change of PE
pdt = 0.5*g*rho0*deedt;

% rate of change of APE
apdt = 0.5*g*rho0*deepdt;

% PE reservoirs
pe_avg = 0.5*g*rho0*eta_avg.*eta_avg;
ape_avg = 0.5*g*rho0*eepa;

% form pdt_check(s)
pdt_check = g*rho0*eta_avg.*dedt;
apdt_check = g*rho0*epa.*depdt;
% and a second version for PE that relies on the divergence
dhu = diff(hudy,1,2);
dhv = diff(hvdx,1,1);
dive = NaN * eta_avg;
dive(2:end-1,2:end-1) = dhu(2:end-1,:) + dhv(:,2:end-1);
dive = dive ./ (G.DX .* G.DY);
pdt_check2 = -g*rho0*eta_avg.*dive;
% dive = divergence of (hu,hv) [m s-1]

% get diagnostics for SW momentum
ub = struct();
vb = struct();
for ii = 1:length(uv_list)
    varname = uv_list{ii};
    uns = ['ubar_',varname];
    ub.(varname) = nc_varget(f_dia,uns);
    ub.(varname)(~G.mask_u) = 0;
    %
    vns = ['vbar_',varname];
    vb.(varname) = nc_varget(f_dia,vns);
    vb.(varname)(~G.mask_v) = 0;
end

% form terms in the energy balance
ue = struct();
ve = struct();
re = struct(); % rho grid
for ii = 1:length(uv_list)
    varname = uv_list{ii};
    ue.(varname) = rho0*hu.*ub.(varname);
    ve.(varname) = rho0*hv.*vb.(varname);
    re.(varname) = ...
        interp2(G.lon_u,G.lat_u,ue.(varname),G.lon_rho,G.lat_rho) + ...
        interp2(G.lon_v,G.lat_v,ve.(varname),G.lon_rho,G.lat_rho);
end

% extra terms to get d(KE)/dt in correct form
kdt_extra_u =  0.5*rho0*hu.*dudt;
kdt_extra_v = 0.5*rho0*hv.*dvdt;
kdt_extra = ...
    interp2(G.lon_u,G.lat_u,kdt_extra_u,G.lon_rho,G.lat_rho) + ...
    interp2(G.lon_v,G.lat_v,kdt_extra_v,G.lon_rho,G.lat_rho) ...
    - 0.5*re.accel;

% check on d(KE)/dt [W m-2] (rho grid)
ke1 = Z_make_SWke(G,rho0,f_his1);
ke2 = Z_make_SWke(G,rho0,f_his2);
kdt_check = (ke2 - ke1)/DT;

% KE reservoir [J m-2] (rho grid)
ke_avg = Z_make_SWke(G,rho0,f_avg);

% output
sw = struct();
% prefix "a" means they apply to the APE budget
%
% fluxes
sw.kdt = re.accel + kdt_extra;
sw.pdt = pdt;
sw.apdt = apdt;
sw.bern = re.xadv + re.yadv + re.cor + re.prsgrd + kdt_extra + sw.pdt;
sw.abern = re.xadv + re.yadv + re.cor + re.prsgrd + kdt_extra + sw.apdt;
sw.sstr = re.sstr;
sw.bstr = re.bstr;
% check fluxes
sw.kdt_check = kdt_check;
sw.pdt_check = pdt_check;
sw.apdt_check = apdt_check;
sw.pdt_check2 = pdt_check2;
% reservoirs
sw.ke = ke_avg;
sw.pe = pe_avg;
sw.ape = ape_avg;

sw_list = {'kdt','pdt','apdt','bern','abern','sstr','bstr', ...
    'kdt_check','pdt_check','apdt_check','pdt_check2','ke','pe','ape'};

for vv = 1:length(sw_list)
    varname = sw_list{vv};
    sw.(varname)(island) = NaN;
end

