% APE_test.m  6/12/2014  Parker MacCready
%
% code to test different methods to calculate APE per unit area
%
% Added lines to make apev positive definite.  There were a few very small
% negative points in the test case (min was -1.8e-9 compared to a median of
% +6 and a mean of +77 in apev).  This could be due to roundoff error.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

dir1 = '/Users/PM5/Documents/roms/output/';

if 0
    basename = 'D2005_his';
    nn = 1836; %
else
    basename = 'salish_2006_4';
    nn = 8701; %
end
dir0 = [dir1,basename,'/'];
%% get hypsometry structure H

ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
fn = [dir0,'ocean_his_',ns,'.nc'];

[G,S,T] = Z_get_basic_info(fn);
island = Z_island(G);
[H] = Z_hyp(G,S,island);
[apea, apea_up, apea_down, zz, rr] = Z_ape(fn,G,S,H);

% also get buoyancy flux
[Fb_a] = Z_make_Fb_a(fn);
Fb10 = Fb_a;
Fb10(Fb10<=0) = 1.e-6;
Fb10 = log10(Fb10);
Fb10 = double(Fb10); %needed for the PNWTOX output
Fb10(island) = nan;

%% plotting
close all
figure
Z_fig(14)
set(gcf,'position',[10 10 1200 1000]);
cax = [2 6];
aa = [-123.25 -122 47 48.4];

subplot(1,2,1)
colormap jet
Z_pcolorcen(G.lon_rho,G.lat_rho,log10(apea));
caxis(cax);
shading flat
colorbar
hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
title(['(a) log_{10}APE_{A} (J m^{-2})'])
Z_dar;
Z_addcoast('combined',Tdir.coast);
axis(aa);
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

subplot(1,2,2)
colormap jet
Z_pcolorcen(G.lon_rho,G.lat_rho,Fb10);
caxis([-4 -1]);
shading flat
colorbar
%hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
title('(b) log_{10}Fb (W m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);
axis(aa);
xlabel('Longitude (deg)')


set(gcf,'PaperPositionMode','auto');
print('-dpng','/Users/PM5/Desktop/PugetSound_APE_Fb.png');





