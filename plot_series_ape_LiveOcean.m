%% Code for plotting APE time series (for LiveOcean).
%

% set up the environment
clear
addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

dir0 = [Tdir.output,'APE_out/'];

load([dir0,'cascadia1_base_lo1_fjord/APE.mat']);

%% plotting
close all
Z_fig(16);
set(gcf,'position',[100 100 1000 650]);
td = ape_ser.td;
ty = 2013 + (td - datenum(2013,1,1))/365;
volume = 5.5531e+12; % volume of the shelf (from info.volume in plot_series.m)
plot(ty,ape_ser.shelf/volume,'-b','linewidth',2);

aa = axis;
axis([2013 2016 aa(3:4)]);

hold on
plot([2014 2014], aa(3:4), '-k')
plot([2015 2015], aa(3:4), '-k')

xlabel('Year')
ylabel('APE (J m^{-3})')
title('SHELF Available Potential Energy');

set(gca,'XTick',[2013.5, 2014.5, 2015.5]);
set(gca,'XTickLabel',{'2013','2014','2015'});

%% printing
if 1
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',[dir0,'APE_131415.png'],'-r0');
end

