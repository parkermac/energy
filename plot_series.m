% Plots the results of make_series with the goal of understanding
% how the energy systems of the shelf and estuary work.

clear
addpath('../alpha/'); Tdir = toolstart;

island_tag = 'salish'; % salish, shelf, abyss, full
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

nfilt = 10; % filter length in days

lw = 3; % linewidth

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
title(['E^{\prime} Rates [W]  ::  ',upper(island_tag),' (',num2str(nfilt),' day filter)'])
legend('dAPE/Dt','dKE/Dt','APE Check','KE check','Bernoulli','Dissipation','Mixing','Background',0)
set(gca,'xlim',dlim);
grid on

%% Choice of APE' Budget or Environmental terms

do_ape = 0;

if do_ape
    subplot(312)
    plot(day,Z_jfilt(P2.pdt',nfilt),'-r', ...
        day,Z_jfilt(P2.rate_check',nfilt),'--k', ...
        day,Z_jfilt(P2.adv',nfilt),'-b', ...
        day,Z_jfilt(P2.diss',nfilt),'-g', ...
        day,Z_jfilt(P2.dpe0dt',nfilt),'-m', ...
        day,Z_jfilt(P2.err',nfilt),'-k', ...
        'linewidth',lw)
    title('APE^{\prime} Rates [W]')
    legend('dAPE/dt','Rate Check','Adv + Conv','Mixing','Background','Error',0)
    set(gca,'xlim',dlim);
    grid on
    
else
    
    % mooring fields in structure M
    load([Tdir.output,'moor_out/cascadia/C',ys,'_RN.mat']);
    M.yd = M.td - day0;
    
    % rivers in structure V
    V = load([Tdir.river,'riverFlow_1998_2009.mat']);
    icr = find(strcmp(V.rname,'fraser'));
    Vmask = V.Qr_year == yn;
    qcr = V.Qr_flow(Vmask,icr);
    
    % Plot Environmental terms
    
    subplot(312)
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
    [xt,yt] = Z_lab('ur');
    text(xt,yt,'Columbia RIVER (3e4 m^{3} s^{-1})','color','b', ...
        'fontweight','bold','fontsize',18, ...
        'horizontalalignment','r')
    
end

%% Reservoirs

subplot(313)
% note the APE terms in P2 do not include the SW part
plot(day,P2.ape,'-r', ...
    day,P2.ape_up,'--r', ...
    day,P2.ape_down,':r', ...
    day,K2.ke - SW.ke,'-b', ...
    day,SW.ke,':c', ...
    day,SW.pe,'-c', ...
    'linewidth',lw)
title('Reservoirs [J]')
legend('APE^{\prime}','APE^{\prime} up','APE^{\prime} down', ...
    'KE^{\prime}','KE^{SW}','PE^{SW}',0)
set(gca,'xlim',dlim);
grid on

xlabel('Days')

%% Printing

if 0
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',[Tdir.output,['energy_out/series_2_',island_tag,'.png']]);
end



