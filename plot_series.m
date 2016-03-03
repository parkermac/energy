% Plots the results of make_series with the goal of understanding
% how the energy systems of the shelf and estuary work.

clear
addpath('../alpha/'); Tdir = toolstart;

island_tag = 'shelf'; % salish, shelf, abyss, full
dlim = [0 365]; % full year is [0 365];
rho0 = 1023.7; % set in the ROMS parameters

pth = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/flux_lp71/';
load([pth,'series_',island_tag,'.mat']); % has island
info.mass = info.volume * rho0;
day = info.td_vec;
ys = datestr(day(1),'yyyy');
yn = str2double(ys);
day0 = datenum(yn,1,1);
day = day - day0;

NT = length(day);

close all
Z_fig(14);
set(gcf,'position',[200 50 1400 900]);
set(gcf,'PaperPositionMode','auto');

nfilt = 30; % filter length in days

lw = 3; % linewidth

if 0
    subplot(311)
    plot(day,Z_jfilt(I2.pdt',nfilt),'-r', ...
        day,Z_jfilt(I2.kdt',nfilt),':r', ...
        day,Z_jfilt(P2.rate_check',nfilt),'--k', ...
        day,Z_jfilt(K2.accel_check'- SW.kdt',nfilt),'.k', ...
        day,Z_jfilt(I2.bern',nfilt),'-b', ...
        day,Z_jfilt(I2.diss' - P2.diss',nfilt),':g', ...
        day,Z_jfilt(P2.diss',nfilt),'-g', ...
        day,Z_jfilt(I2.background',nfilt),'-m', ...
        'linewidth',lw)
    title(['(a) E^{\prime} Rates [W]  ::  ',upper(island_tag),' (',num2str(nfilt),' day filter)'])
    legend('dAPE/Dt','dKE/Dt','APE Check','KE check','Bernoulli','Dissipation','Mixing','Background',0)
    set(gca,'xlim',dlim);
    grid on
        
else
    subplot(311)
    plot(day,Z_jfilt(P2.pdt',nfilt)/info.mass,'--b', ...
        day,Z_jfilt(P2.adv',nfilt)/info.mass,'-r', ...
        day,Z_jfilt(P2.diss',nfilt)/info.mass,'-m', ...
        day,Z_jfilt(P2.dpe0dt',nfilt)/info.mass,'-g', ...
        day,Z_jfilt(P2.pdt' - P2.rate_check',nfilt)/info.mass,'-k', ...
        'linewidth',lw)
    title(['(a) APE^{\prime} Rates [W kg^{-1}]  ::  ',upper(island_tag),' (',num2str(nfilt),' day filter)'])
    legend('dAPE/dt','Adv + Conv','Mixing','Background','Error',0)
    set(gca,'xlim',dlim);
    grid on
    
end

%% Plot Environmental terms

% mooring fields in structure M

load([Tdir.output,'moor_out/cascadia/C',ys,'_RN.mat']);
M.yd = M.td - day0;

% rivers in structure V
if strcmp(island_tag,'shelf')
    rivername = 'Columbia';
    rivscale = 6e4;
else
    rivername = 'Fraser'
    rivscale = 3e4;
end
V = load([Tdir.river,'riverFlow_1998_2009.mat']);
icr = find(strcmp(V.rname,lower(rivername)));
Vmask = V.Qr_year == yn;
qcr = V.Qr_flow(Vmask,icr);

subplot(312)
hold on
if 1
    wind = Z_8d(M.svstr')';
    wday = M.yd;
else
    wind = Z_godin(M.svstr')';
    wind = wind(36:24:8760-36);
    wind = Z_jfilt(wind',nfilt)';
    wday = day;
end
wind(isnan(wind)) = 0;
mask = wind > 0;
windp = wind; windp(~mask) = 0;
windm = wind; windm(mask) = 0;
fill([dlim(1) wday dlim(2)],[0 windp 0],'g','edgecolor','g')
fill([dlim(1) wday dlim(2)],[0 windm 0],'m','edgecolor','m')
title('(b) Wind and River'); xlim(dlim); grid on
ylim([-.1 .3]);
[xt,yt] = Z_lab('ul');
text(xt,yt,'DOWNWELLING Wind Stress [Pa]','color','g','fontweight','bold','fontsize',18)
[xt,yt] = Z_lab('ll');
text(xt,yt,'UPWELLING Wind Stress [Pa]','color','m','fontweight','bold','fontsize',18)

plot([1:365],qcr/rivscale,'-b','linewidth',3)
[xt,yt] = Z_lab('ur');
text(xt,yt,[rivername,' River [',num2str(rivscale),' m^{3} s^{-1}]'],'color','b', ...
    'fontweight','bold','fontsize',18, ...
    'horizontalalignment','r')


%% Reservoirs

subplot(313)
% note the APE terms in P2 do not include the SW part
plot(day,P2.ape/info.mass,'-r', ...
    day,P2.ape_up/info.mass,'--r', ...
    day,P2.ape_down/info.mass,':r', ...
    day,(K2.ke - SW.ke)/info.mass,'-b', ...
    day,SW.ke/info.mass,':c', ...
    day,SW.pe/info.mass,'-c', ...
    'linewidth',lw)
title('(c) Reservoirs [J kg^{-1}]')
legend('APE^{\prime}','APE^{\prime} up','APE^{\prime} down', ...
    'KE^{\prime}','KE^{SW}','PE^{SW}',0)
set(gca,'xlim',dlim);
grid on

xlabel('Days')

%% Printing

if 1
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',[Tdir.output,['energy_out/series_new_',island_tag,'.png']]);
end

%% make a second figure

w = wind(36:24:8760 - 36);

figure;
plot(w,Z_jfilt(P2.adv',nfilt)/info.mass,'*r');
grid on
xlabel('8-day Filtered Meridional Wind Stress [Pa]')
ylabel('APE^{\prime} Rate: Adv+Conv [W kg^{-1}]')
axis square
aa = axis;
hold on
plot(aa(1:2),[0,0],'-k')
plot([0,0],aa(3:4),'-k')
axis(aa)

%% Printing

if 1
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',[Tdir.output,['energy_out/adv_conv_vs_wind_',island_tag,'.png']]);
end



