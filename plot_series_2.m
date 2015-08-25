% plot_series_2.m  8/21/2015  Parker MacCready
%
% plots the results of make_series with the goal of understanding
% how the energy systems of the shelf and estuary work

clear
addpath('../alpha/'); Tdir = toolstart;

island_tag = 'shelf'; % salish, shelf, abyss, full
dlim = [0 365]; % full year is [0 365];

pth = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/flux_lp71/';
load([pth,'series_',island_tag,'.mat']); % has island
day = info.td_vec;
ys = datestr(day(1),'yyyy');
yn = str2double(ys);
day0 = datenum(yn,1,1);
day = day - day0;

NT = length(day);

close all
Z_fig(14);
set(gcf,'position',[200 50 1800 1200]);
set(gcf,'PaperPositionMode','auto');

nfilt = 30;

subplot(311)
plot(day,Z_jfilt(I2.edt',nfilt),'-r', ...
    day,Z_jfilt(I2.kdt',nfilt),'--r', ...
    day,Z_jfilt(I2.pdt',nfilt),':r', ...
    day,Z_jfilt(I2.bern',nfilt),'-b', ...
    day,Z_jfilt(I2.diss',nfilt),'-g', ...
    day,Z_jfilt(K2.wind_work',nfilt),':g', ...
    day,Z_jfilt(I2.err',nfilt),'-k', ...
    'linewidth',2)
title(['Rates (E^{\prime}) [W] ',upper(island_tag)])
legend('edt','kdt','pdt','bern','diss','background','err',0)
set(gca,'xlim',dlim);
grid on

subplot(312)
plot(day,Z_jfilt(P2.pdt',nfilt),'-r', ...
    day,Z_jfilt(P2.adv',nfilt),'-b', ...
    day,Z_jfilt(P2.diss',nfilt),'-g', ...
    day,Z_jfilt(P2.dpe0dt',nfilt),'-m', ...
    day,Z_jfilt(P2.err',nfilt),'-k', ...
    'linewidth',2)
title(['Rates (APE^{\prime}) [W] ',upper(island_tag)])
legend('pdt','adv','diff','background','err',0)
set(gca,'xlim',dlim);
grid on

subplot(313)
% note I think the APE terms do not include the SW part
plot(day,P2.ape,'-r', ...
    day,P2.ape_up,'--r', ...
    day,P2.ape_down,':r', ...
    day,K2.ke - SW.ke,'-b', ...
    day,SW.ke,':c', ...
    day,SW.pe,'-c', ...
    'linewidth',2)
title('Reservoirs [J]')
legend('APE^{\prime}','APE^{\prime} up','APE^{\prime} down', ...
    'KE^{\prime}','KE^{SW}','PE^{SW}',0)
set(gca,'xlim',dlim);
grid on









