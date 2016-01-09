% Plots the APEv field with sections to show depth structure.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');
fn = '/Users/PM5/Documents/roms/output/pnwtox_2006/ocean_his_5000.nc';

%% call the function
[G,S,T] = Z_get_basic_info(fn);
island = Z_island(G);
[H] = Z_hyp(G,S,island);

%% calculate APE

%% calculate APE

% constants
g = 9.81;
offset = 1000; % to add to the density anomaly

% space info
eta = nc_varget(fn,'zeta');
[zr,zw] = Z_s2z(G.h,0*eta,S);
DZ = diff(zw);

salt = nc_varget(fn,'salt');
temp = nc_varget(fn,'temp');

[~,den1,~,~] = Z_roms_eos_vec(salt,temp,zr);
rho = den1 + offset; % use potential density

[D] = Z_flat(rho,rho,zr,zw,H);

zz = D.Z - D.Zf;

% create exact APE (neglecting compressibility)
F = g*(D.intRf - D.intRfzf);
apev = g*zz.*rho - F;
apev(apev <= 0) = 1e-10;
apea = squeeze(sum(apev.*DZ)); % [J m-2]

%% plotting

close all
figure
Z_fig(14)
set(gcf,'position',[10 10 1300 1000]);

%% map

lat_list = [49.2, 47, 44];

ax_mat = [-124.2 -122.8 -400 0; -125.5 -124 -400 0; -125.5 -124 -400 0];

subplot(121);
colormap jet
Z_pcolorcen(G.lon_rho,G.lat_rho,log10(apea));
caxis([2 6]);
shading flat
colorbar('south')
hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
title('log10 APE_{A} (J m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')


for ii = 1:length(lat_list)
    lat = lat_list(ii);
    jlat = find(G.lat_rho(:,1) >= lat); jlat = jlat(1);
    llat = G.lat_rho(jlat,1);
    hold on
    plot([ax_mat(ii,1) ax_mat(ii,2)],[llat llat], '-k','linewidth',3)
end

%% sections


for ii = 1:length(lat_list)
    lat = lat_list(ii);
    jlat = find(G.lat_rho(:,1) >= lat); jlat = jlat(1);
    Lat = G.lat_rho(jlat,:);
    Lon = G.lon_rho(jlat,:);
    Apev = squeeze(apev(:,jlat,:));
    Sig = squeeze(den1(:,jlat,:));
    Zr = squeeze(zr(:,jlat,:));
    Zw = squeeze(zw(:,jlat,:));
    
    Zf = [Zw(1,:) ; Zr ; Zw(end,:)];
    Apevf = [Apev(1,:) ; Apev ; Apev(end,:)];
    Sigf = [Sig(1,:) ; Sig ; Sig(end,:)];
    Lonf = ones(size(Zf,1),1) * Lon;
    
    subplot(length(lat_list),2,2*ii)
    pcolor(Lonf,Zf,log10(Apevf))
    caxis([0 4])
    shading interp
    colorbar('east')
    
    hold on
    contour(Lonf, Zf, Sigf,[20:.5:30],'-k')
    
    plot(Lon,Zw(1,:),'-k','linewidth',3)
    
    axis(ax_mat(ii,:))
    
    xlabel('Longitude (deg)')
    ylabel('Z (m)')
    if ii == 1
        title('log10 APE_{V} (J m^{-3})')
    end
    
    box on
end


