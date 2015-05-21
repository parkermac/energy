%% flux_test.m  6/12/2014  Parker MacCready
%
% does area integrations of the vertically-integrated flux terms
% calculated by flux.m (via flux_driver.m)

clear; addpath('../alpha/'); Tdir = toolstart;
odir_top = [Tdir.output,'energy_out/'];
addpath('./Zfun');
close all

% set the file name to look at
nn = 1836;
ncpad = ['0000',num2str(nn)]; ncpad = ncpad(end-4:end);
infile = ['flux_',ncpad,'.mat'];

% what to look at
do_flux_plot = 1;
do_res_plot = 1;
show_balance = 1;
check_mask = 1;
do_lp = 1;

[basename,nn_vec,dir0] = Z_runspec_raw;

if do_lp
    idir = [odir_top,basename,'/flux_lp71/'];
else
    idir = [odir_top,basename,'/flux_raw/'];
end

load([idir,'G.mat']); % has island
island_orig = G.island;
ISN = sum(island_orig(:)); % note how many land points island has

%% load file

inname = [idir,infile];
load(inname);

% NOTE 6/9/2014 There are a few (21) extra NaN's in the KE budget associated
% with river mouths.  For integration purposes we will get rid of these
% using the following two lines
G.island = isnan(k2.accel_check); % redefine G.island
DA = G.DX .* G.DY;
DA(G.island) = NaN;

%% get derived quantities
[p2,k2,e2,sw,info] = Z_make_derived(p2,k2,sw,info);

%% check the mask
if check_mask
    Z_check_mask(G,info,p2,k2,sw,ISN,island_orig,Tdir);
end

%% show balances (area integrals)
if show_balance
    Z_show_balance(p2,k2,e2,sw,info,DA);
end

%%  flux plot
if do_flux_plot
    Z_flux_plot(G,p2,k2,e2,sw,do_lp,Tdir);
end

%% reservoir plotting
if do_res_plot
    Z_res_plot(G,p2,k2,e2,sw,do_lp,Tdir,info);
end







