% plot_series.m  6/18/2014  Parker MacCready
%
% plots the results of make_series

clear
addpath('../alpha/'); Tdir = toolstart;

island_tag = 'full'; % salish, shelf, abyss, full
dlim = [0 365]; % full year is [0 365];

pth = '/Users/PM3/Documents/tools_output/energy_out/Cdia2005/flux_raw/';
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


% low pass filtering
if 1
    for vv = 1:length(info.ts_list)
        vname = info.ts_list{vv};
        P2.(vname) = Z_godin(P2.(vname)')';
    end
    for vv = 1:length(info.uv_list)
        vname = info.uv_list{vv};
        K2.(vname) = Z_godin(K2.(vname)')';
    end
    for vv = 1:length(info.e2_list)
        vname = info.e2_list{vv};
        E2.(vname) = Z_godin(E2.(vname)')';
    end
    for vv = 1:length(info.sw_list)
        vname = info.sw_list{vv};
        SW.(vname) = Z_godin(SW.(vname)')';
    end  
end


subplot(411)
plot(day,E2.edt - SW.edt,'-r', ...
    day,E2.bern - SW.bern,'-g', ...
    day,E2.diss - SW.diss,'-b', ...
    day,E2.err - SW.err,'-k', ...
    'linewidth',2)
title(['Rates (Full-SW) [W] ',upper(island_tag)])
legend('dE/dt','Bernoulli','Diss','Err',0)
set(gca,'xlim',dlim);
grid on

subplot(412)
plot(day,P2.hdiff,'-r', ...
    day,P2.vdiff,'-b', ...
    day,SW.sstr,'-g', ...
    'linewidth',2)
title('Selected Turbulent Rates [W]')
legend('PE hdiff','PE vdiff','SW sstr',0)
set(gca,'xlim',dlim);
grid on

subplot(413)
plot(day,P2.ape,'-r', ...
    day,K2.ke - SW.ke,'-b', ...
    day,SW.ape,'-m', ...
    day,SW.ke,'-c', ...
    'linewidth',2)
title('Reservoirs [J]')
legend('APE','KE (Full-SW)','APE SW','KE SW',0)
set(gca,'xlim',dlim);
grid on

% % Separate PW and KE rates, and checks on them
% subplot(325)
% plot(day,P2.rate - SW.apdt,'-r', ...
%     day,P2.rate_check - SW.apdt_check,'-m', ...
%     day,K2.accel - SW.kdt,'-b', ...
%     day,K2.accel_check - SW.kdt_check,'-c', ...
%     'linewidth',2)
% title('Rates (Full-SW) [W]')
% legend('APE','APE check','KE','KE check',0)
% set(gca,'xlim',dlim);
% grid on

%% mooring fields in structure M
load([Tdir.output,'moor_out/C',ys,'_RN.mat']);
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







