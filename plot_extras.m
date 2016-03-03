%% plotting extras

indir = '/Users/PM5/Documents/tools_output/energy_out/D2005/extras_lp71/';

load([indir,'extras_01836.mat']);
load([indir,'G.mat']);

close all
figure
Z_fig(16)
set(gcf,'position',[10 10 1300 1000]);
colormap jet

Z_pcolorcen(G.lon_rho,G.lat_rho, ex2.wind_work);
shading flat
colorbar('eastoutside')
title('Wind work (W m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);

[xt,yt] = Z_lab('lr');
text(xt,yt,datestr(info.td_avg,'dd-mmm-YYYY'), ...
    'horizontalalignment','r');


xlabel('Longitude (deg)')
ylabel('Latitude (deg)')


