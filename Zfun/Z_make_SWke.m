function ke = Z_make_SWke(G,rho0,fn)
% 6/12/2014 Parker MacCready
%
% returns vertically-integrated SW KE [J m-2]
u = nc_varget(fn,'ubar');
v = nc_varget(fn,'vbar');
eta = nc_varget(fn,'zeta');

uu = u.*u;
vv = v.*v;

uu(~G.mask_u) = 0;
vv(~G.mask_v) = 0;

uuvv = interp2(G.lon_u,G.lat_u,uu,G.lon_rho,G.lat_rho) + ...
        interp2(G.lon_v,G.lat_v,vv,G.lon_rho,G.lat_rho);
    
ke = 0.5*rho0*(G.h + eta).*uuvv;