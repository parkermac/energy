% APE_func_test.m  6/12/2014  Parker MacCready
%
% code to test Z_ape

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');
fn = '/Users/PM5/Documents/roms/output/D2005_avg/ocean_avg_1836.nc';

%% call the function
[G, apea, apea_up, apea_down, zz, rr] = Z_ape(fn);

%% plotting
close all
figure
Z_fig(14)
set(gcf,'position',[10 10 2500 1000]);
cax = [2 6];
fld_list = {'apea','apea_up','apea_down'};
for ii = 1:length(fld_list)
    fld_name = fld_list{ii};
    eval(['fld = ',fld_name,';']);
    subplot(1,3,ii)
    colormap jet
    Z_pcolorcen(G.lon_rho,G.lat_rho,log10(fld));
    caxis(cax);
    shading flat
    colorbar
    hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
    title(['log10 ',strrep(fld_name,'_',' '),' (J m^{-2})'])
    Z_dar;
    Z_addcoast('combined',Tdir.coast);
end


