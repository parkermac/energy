% plot_series.m  6/18/2014  Parker MacCready
%
% plots the results of make_series

clear
addpath('../alpha/'); Tdir = toolstart;

island_tag = 'full'; % salish, shelf, abyss, full
dlim = [0 365]; % full year is [0 365];

pth = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/flux_lp71/';
load([pth,'series_',island_tag,'.mat']); % has island
day = info.td_vec;
ys = datestr(day(1),'yyyy');
yn = str2num(ys);
day0 = datenum(yn,1,1);
day = day - day0;

NT = length(day);

close all
Z_fig(14);
set(gcf,'position',[200 50 1800 1200]);
set(gcf,'PaperPositionMode','auto');

subplot(411)
plot(day,I2.edt,'-r', ...
    day,I2.pdt,'--r', ...
    day,I2.kdt,':r', ...
    day,I2.background,'-m', ...
    day,I2.bern,'-g', ...
    day,I2.diss,'-b', ...
    day,I2.err,'-k', ...
    'linewidth',2)
title(['Rates (Full-SW) [W] ',upper(island_tag)])
legend('dE/dt','pdt','kdt','background','Bernoulli','Diss','Err',0)
set(gca,'xlim',dlim);
grid on

subplot(412)
plot(day,P2.hdiff,'-r', ...
    day,P2.vdiff,'-b', ...
    day,SW.sstr,'-g', ...
    day,K2.wind_work,'m', ...
    'linewidth',2)
title('Selected Turbulent Rates [W]')
legend('PE hdiff','PE vdiff','SW sstr','Wind Work',0)
set(gca,'xlim',dlim);
grid on

subplot(413)
plot(day,P2.ape,'-r', ...
    day,P2.ape_up,'--r', ...
    day,P2.ape_down,':r', ...
    day,K2.ke - SW.ke,'-b', ...
    day,SW.pe,'-m', ...
    day,SW.ke,'-c', ...
    'linewidth',2)
title('Reservoirs [J]')
legend('APE','APE up','APE down','KE (Full-SW)','APE SW','KE SW',0)
set(gca,'xlim',dlim);
grid on

%% mooring fields in structure M
load([Tdir.output,'moor_out/cascadia/C',ys,'_RN.mat']);
M.yd = M.td - day0;

%% rivers in structure V
V = load([Tdir.river,'riverFlow_1998_2009.mat']);
icr = find(strcmp(V.rname,'fraser'));
Vmask = V.Qr_year == yn;
qcr = V.Qr_flow(Vmask,icr);

%%

subplot(414)
%plot(M.yd,Z_jfilt(M.svstr',40),'-k')
hold on
w8 = Z_8d(M.svstr')';
w8(isnan(w8)) = 0;
mask = w8 > 0;
w8p = w8; w8p(~mask) = 0;
w8m = w8; w8m(mask) = 0;
%plot(M.yd,w8,'-k')
fill([dlim(1) M.yd dlim(2)],[0 w8p 0],'g','edgecolor','g')
fill([dlim(1) M.yd dlim(2)],[0 w8m 0],'m','edgecolor','m')
title('Low-Passed NS Windstress [Pa]'); xlim(dlim); grid on
ylim([-.1 .3]);
[xt,yt] = Z_lab('ul');
text(xt,yt,'DOWNWELLING','color','g','fontweight','bold','fontsize',18)
[xt,yt] = Z_lab('ll');
text(xt,yt,'UPWELLING','color','m','fontweight','bold','fontsize',18)

plot([1:365],qcr/3e4,'-b','linewidth',3)

%%
print('-djpeg100',['/Users/PM5/Desktop/',island_tag,'.jpg']);







