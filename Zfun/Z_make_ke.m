function ke = Z_make_ke(G,S,rho0,fn)
% 6/12/2014 Parker MacCready
% creates vertically-integrated KE on the rho grid for a given filename

% first create the cell thickness terms on the u- and v-grids
h_u = interp2(G.lon_rho,G.lat_rho,G.h,G.lon_u,G.lat_u);
h_v = interp2(G.lon_rho,G.lat_rho,G.h,G.lon_v,G.lat_v);

eta = nc_varget(fn,'zeta');

% note the unusual interpolation needed to avoid introducing
% too many NaN's
eta_u = NaN * G.lon_u;
for ii = 1:size(G.lon_rho,1)
    eta_u(ii,:) = interp1(G.lon_rho(ii,:),eta(ii,:),G.lon_u(ii,:));
end
eta_v = NaN * G.lon_v;
for ii = 1:size(G.lon_rho,2)
    eta_v(:,ii) = interp1(G.lat_rho(:,ii),eta(:,ii),G.lat_v(:,ii));
end

zw_u = Z_s2z_w(h_u,eta_u,S); dz_u = diff(zw_u);
zw_v = Z_s2z_w(h_v,eta_v,S); dz_v = diff(zw_v);

clear zw_u zw_v

u = nc_varget(fn,'u');
v = nc_varget(fn,'v');

ke3_u = 0.5 * rho0 * u .* u .* dz_u;
ke3_v = 0.5 * rho0 * v .* v .* dz_v;

clear u v dz_u dz_v
clear h_u h_v

ke_u = squeeze(sum(ke3_u));
ke_v = squeeze(sum(ke3_v));

clear ke3_u ke3_v

ke_u(~G.mask_u) = 0;
ke_v(~G.mask_v) = 0;

ke = ...
    interp2(G.lon_u,G.lat_u,ke_u,G.lon_rho,G.lat_rho) + ...
    interp2(G.lon_v,G.lat_v,ke_v,G.lon_rho,G.lat_rho);