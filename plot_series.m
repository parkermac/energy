%% Plots the results of make_series with the goal of understanding
% how the energy systems of the shelf and estuary work.

% Brings in other data such as Wind Work, TEF, Rivers, and Wind.

clear
addpath('../alpha/'); Tdir = toolstart;

% choices
island_tag = 'salish'; % salish, shelf, abyss, full
nfilt = 5; % filter length in days
yh_line = 4800; % hour for a marker (like a line)
yd_line = yh_line/24; % yearday for the same marker
do_print = 0;

% set things that are unlikely to change
dlim = [0 365]; % full year is [0 365];
rho0 = 1023.7; % set in the ROMS parameters

%% Load the energy budget results, with structures of 2D terms:
%
% P2, K2, SW: Basic energy budget terms
% E2: Derived complete energy budget terms
% I2: Derived internal energy budget terms
% EX2: Extra terms (geostrophic wind work)
%
% each 363 daily times
% starting 02-Jan-2005 11:30:00, ending 30-Dec-2005 11:30:00
% presumably the missing days are due to the requirements of the Godin
% filter.

pth = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/flux_lp71/';
load([pth,'series_',island_tag,'.mat']); % has island
info.mass = info.volume * rho0; % approximate mass [kg] of the volume
day = info.td_vec; % datenum time vector
ys = datestr(day(1),'yyyy');
yn = str2double(ys);
day0 = datenum(yn,1,1);
day = day - day0; % now this is yearday
NT = length(day);

%% Load TEF terms, from Sarah Giddings
%
% 8761 hourly fields
% TEF.tdQ is yearday

TEF = load([Tdir.output, ...
    'energy_out/Cdia2005_TEF/Cdia2005_his_his_EH1_SimpleSection.mat']);

%% Load mooring fields in structure M
%
% 8760 hourly fields

load([Tdir.output,'moor_out/cascadia/C',ys,'_RN.mat']);
M.yd = M.td - day0;

M.wind = Z_8d(M.svstr')';
% note Z_8d must only be used with hourly data
M.wday = M.yd; % hourly yearday vector
wind = M.wind;
wind(isnan(wind)) = 0;
mask = wind > 0;
M.windp = wind; M.windp(~mask) = 0;
M.windm = wind; M.windm(mask) = 0;

%% Load rivers in structure R
%
% 365 daily values

R = load([Tdir.river,'riverFlow_1998_2009.mat']);
Rmask = R.Qr_year == yn;
iriv = find(strcmp(R.rname,lower('Columbia')));
R.qr_columbia = R.Qr_flow(Rmask,iriv);
iriv = find(strcmp(R.rname,lower('Fraser')));
R.qr_fraser = R.Qr_flow(Rmask,iriv);
R.yd = 1:365; % a synthetic yearday vector

%% PLOTTING

close all
Z_fig(18);

%% Plot Environmental terms

set(gcf,'position',[50 50 1400 900]);
set(gcf,'PaperPositionMode','auto');
fs = 18; % fontsize

hold on

subplot(311)
fill([dlim(1) M.wday dlim(2)],[0 M.windp 0],'r','edgecolor','r')
hold on
fill([dlim(1) M.wday dlim(2)],[0 M.windm 0],'b','edgecolor','b')
xlim(dlim)
ylim([-.1 .3]);
[xt,yt] = Z_lab('ul');
text(xt,yt,'(a) Wind Stress [Pa] 8-Day Filter','fontweight','bold','fontsize',fs)
[xt,yt] = Z_lab('ur');
text(xt,yt-.05,'DOWNWELLING','color','r','fontweight','bold','fontsize',fs,'horizontalalignment','r')
[xt,yt] = Z_lab('lr');
text(xt,yt,'UPWELLING','color','b','fontweight','bold','fontsize',fs,'horizontalalignment','r')
% add a line at a specific time
aa = axis; hold on; plot(yd_line*[1 1],aa(3:4),'-k'); axis(aa);

subplot(312)
plot(R.yd,R.qr_columbia/1000,'-g','linewidth',3)
hold on
plot(R.yd,R.qr_fraser/1000,'-m','linewidth',3)
xlim(dlim)
ylim([0 12])
[xt,yt] = Z_lab('ul');
text(xt,yt,'(b) Daily River Flow [1000 m^{3} s^{-1}]','fontweight','bold','fontsize',fs)
[xt,yt] = Z_lab('ur');
text(xt,yt,'Columbia River','color','g', ...
    'fontweight','bold','fontsize',fs,'horizontalalignment','r')
text(xt,yt-2,'Fraser River','color','m', ...
    'fontweight','bold','fontsize',fs,'horizontalalignment','r')
% add a line at a specific time
aa = axis; hold on; plot(yd_line*[1 1],aa(3:4),'-k'); axis(aa);

subplot(313)
plot(TEF.tdQ,Z_jfilt(TEF.Qin/1000,24*nfilt),'-k','linewidth',3)
xlim(dlim)
ylim([0 300])
[xt,yt] = Z_lab('ul');
text(xt,yt,'(c) Juan de Fuca Total Exchange Flow [1000 m^{3} s^{-1}]', ...
    'fontweight','bold','fontsize',fs)
% add a line at a specific time
aa = axis; hold on; plot(yd_line*[1 1],aa(3:4),'-k'); axis(aa);
xlabel('Days')

% Printing
if do_print
    print('-dpng',[Tdir.output,'energy_out/env_series_',island_tag,'.png']);
end


%% Plot Energy Terms

figure;

set(gcf,'position',[150 150 1400 900]);
set(gcf,'PaperPositionMode','auto');
fs = 18;
lw = 3; % linewidth

% Reservoirs

subplot(211)
% note the APE terms in P2 do not include the SW part
plot(day,P2.ape/info.mass,'-r', ...
    day,P2.ape_up/info.mass,'--r', ...
    day,P2.ape_down/info.mass,':r', ...
    day,(K2.ke - SW.ke)/info.mass,'-b', ...
    day,SW.ke/info.mass,':c', ...
    day,SW.pe/info.mass,'-c', ...
    'linewidth',lw)
xlim(dlim)
[xt,yt] = Z_lab('ul');
text(xt,yt,['(a) Reservoirs [J kg^{-1}] : ',upper(island_tag)], ...
    'fontweight','bold','fontsize',fs)
legend('APE^{\prime}','APE^{\prime} up','APE^{\prime} down', ...
    'KE^{\prime}','KE^{SW}','PE^{SW}','location','east')
% add a line at a specific time
aa = axis; hold on; plot(yd_line*[1 1],aa(3:4),'-k'); axis(aa);

% budget terms

subplot(212)

scl = 1e6;
plot(day,scl*Z_jfilt(P2.pdt',nfilt)/info.mass,'-r', ...
    day,scl*Z_jfilt(P2.adv',nfilt)/info.mass,'-b', ...
    day,scl*Z_jfilt(P2.diss',nfilt)/info.mass,'-m', ...
    day,scl*Z_jfilt(P2.dpe0dt',nfilt)/info.mass,'-g', ...
    day,scl*Z_jfilt(P2.pdt' - P2.rate_check',nfilt)/info.mass,'-k', ...
    'linewidth',lw)
legend('dAPE/dt','Adv + Conv','Mixing','Background','Error','location','northeast')
[xt,yt] = Z_lab('ul');
text(xt,yt,['(b) APE^{\prime} Rates [10^{-6} W kg^{-1}]', ...
    ' ',num2str(nfilt),'-Day Filter'], ...
    'fontweight','bold','fontsize',fs)
xlim(dlim)
% add a line at a specific time
aa = axis; hold on; plot(yd_line*[1 1],aa(3:4),'-k'); axis(aa);
xlabel('Days')

% Printing
if do_print
    print('-dpng',[Tdir.output,'energy_out/en_series_',island_tag,'.png']);
end

%% Corrolations

figure;
set(gcf,'position',[300 300 1400 600]);
set(gcf,'PaperPositionMode','auto');

subplot(1,3,[1 2])

scl = 1e6;
plot(day,scl*Z_jfilt(EX2.wind_work',nfilt)/info.mass,':b', ...
    day,scl*Z_jfilt(P2.adv',nfilt)/info.mass,'-b', ...
    'linewidth',lw)
legend('Geostrophic Wind Work','Adv + Conv','location','northeast')
[xt,yt] = Z_lab('ul');
text(xt,yt,'(a)','fontweight','bold','fontsize',fs)
xlim(dlim)
grid on
xlabel('Days')
ylabel('[10^{-6} W kg^{-1}]')

subplot(133)

w = M.wind(36:24:8760 - 36);
mask_upwelling = w < 0;
mask_downwelling = w >= 0;

scl = 1e6;
lh1 = plot(scl*Z_jfilt(EX2.wind_work(mask_upwelling)',nfilt)/info.mass, ...
    scl*Z_jfilt(P2.adv(mask_upwelling)',nfilt)/info.mass,'*b');
hold on
lh2 = plot(scl*Z_jfilt(EX2.wind_work(mask_downwelling)',nfilt)/info.mass, ...
    scl*Z_jfilt(P2.adv(mask_downwelling)',nfilt)/info.mass,'*r');
lh3 = plot([0 .3],[0 .15],'-k');
axis([-.05 .3 -.05 .15]);
aa = axis;
[xt,yt] = Z_lab('ul');
text(xt,yt,'(b)','fontweight','bold','fontsize',fs)
grid on
legend([lh1;lh2;lh3],'Upwelling','Downwelling','1:2 Line','location','north')
xlabel('Geostrophic Wind Work [10^{-6} W kg^{-1}]')
ylabel('APE^{\prime} Rate: Adv+Conv [10^{-6} W kg^{-1}]')
hold on
plot(aa(1:2),[0,0],'-k')
plot([0,0],aa(3:4),'-k')
axis(aa)

% Printing
if do_print
    print('-dpng',[Tdir.output,['energy_out/adv_conv_vs_wind_',island_tag,'.png']]);
end


%% More correlations

figure;
set(gcf,'position',[300 300 1400 600]);
set(gcf,'PaperPositionMode','auto');

scl1 = 0.5;
scl2 = 1e6;
plot(day,scl1*Z_jfilt(SW.ke',nfilt)/info.mass,'r', ...
    day,scl2*Z_jfilt(P2.hdiff' + P2.vdiff',nfilt)/info.mass,'-b', ...
    'linewidth',lw)
legend('1/2 SW KE [J kg^{-1}]','Mixing [10^{-6} W kg^{-1}]','location','northeast')
xlim(dlim)
grid on
xlabel('Days')
ylabel('')

% Printing
if do_print
    print('-dpng',[Tdir.output,['energy_out/newcorr_',island_tag,'.png']]);
end




%% extra

% tef_qin = Z_jfilt(TEF.Qin,nfilt*24);
% tef_day = TEF.tdQ;
% 
% tef_qin = tef_qin(36:24:8760 - 36);
% tef_day = tef_day(36:24:8760 - 36);
% 
% figure;
% plot(tef_qin, ...
%     Z_jfilt(P2.adv',nfilt)/info.mass,'*r');
% grid on
% xlabel('Qin [m^{3} s^{-1}]')
% ylabel('APE^{\prime} Rate: Adv+Conv [W kg^{-1}]')
% aa = axis;
% hold on
% plot(aa(1:2),[0,0],'-k')
% plot([0,0],aa(3:4),'-k')
% axis(aa)
% 
% % Printing
% set(gcf,'position',[50 50 600 600]);
% set(gcf,'PaperPositionMode','auto');
% print('-dpng',[Tdir.output,['energy_out/adv_conv_vs_qin_',island_tag,'.png']]);
% 






