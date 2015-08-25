%% Code for plotting APE time series (for SciDAC).
%

% set up the environment
clear
addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

dir0 = [Tdir.output,'APE_out/'];

bn_list = {'S45_OcAcT2005R2005_2005','S45_OcAcT2005R2005_2099', ...
    'S85_OcAcT2005R2005_2005','S85_OcAcT2005R2005_2099'};

bb_list = {'S45_2005','S45_2099','S85_2005','S85_2099'};

for ii = 1:length(bn_list)
    bn = bn_list{ii};
    bb = bb_list{ii};
    a.(bb) = load([dir0,bn,'/APE.mat']);
end

b = load([Tdir.output,'energy_out/Cdia2005/flux_lp71/series_shelf.mat']);

%% plotting
close all
Z_fig(16);
set(gcf,'position',[100 100 1000 650]);
c_list = {'b','r','b','r'};
ls_list = {'--','--','-','-'};
bb_list_2 = {};
for ii = 1:length(bb_list)
    bb = bb_list{ii};
    bb_list_2{ii} = strrep(bb,'_',' ');
    td = a.(bb).ape_ser.td;
    td = td - td(1); % days from the start of the year
    ape = a.(bb).ape_ser.shelf;
    lh(ii) = plot(td,ape,[c_list{ii},ls_list{ii}],'linewidth',2);
    hold on
end
% now add the Cascadia 2005 run:
td = b.info.td_vec;
td = td - td(1);
ape = b.P2.ape;
lh(5) = plot(td,ape,'-g','linewidth',2);
bb_list_2{5} = 'Cascadia 2005';

aa = axis;
axis([0 365 aa(3:4)]);

legend(lh',bb_list_2,0);
xlabel('Days')
ylabel('APE (J)')
title('SHELF Available Potential Energy');

set(gcf,'PaperPositionMode','auto');
print('-dpng',[dir0,'APE_plot_1.png']);

